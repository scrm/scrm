/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013, 2014 Paul R. Staab, Sha (Joe) Zhu, Dirk Metzler and Gerton Lunter
 * 
 * This file is part of scrm.
 * 
 * scrm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "time_interval.h"

/* --------------------------------------------------------------------
 * TimeInterval
 * -------------------------------------------------------------------*/

TimeInterval::TimeInterval() {
  this->tii_ = NULL;
  this->start_height_ = 0;
  this->end_height_ = 0;
}


TimeInterval::TimeInterval(TimeIntervalIterator* tii, double start_height, double end_height){
  assert( tii != NULL );
  this->tii_ = tii;
  this->start_height_ = start_height;
  this->end_height_ = end_height;
}


/* --------------------------------------------------------------------
 * TimeIntervalIterator
 * -------------------------------------------------------------------*/

TimeIntervalIterator::TimeIntervalIterator(Forest *forest) {
  // Used only for unit testing, and hence private.
  this->forest_ = forest;
  this->contemporaries_ = &(forest->contemporaries_);
  this->contemporaries_->clear();
  this->node_iterator_ = forest->nodes()->iterator();
  this->good_ = false;
  this->inside_node_ = NULL;
  this->current_time_ = 0; 
  forest->writable_model()->resetTime();
}

TimeIntervalIterator::TimeIntervalIterator(Forest* forest, 
                                           Node* start_node) {
  this->forest_ = forest;
  this->contemporaries_ = &forest->contemporaries_;
  this->model_ = forest->model_;

  this->good_ = true;
  this->inside_node_ = NULL;
  this->node_iterator_ = forest->nodes()->iterator(start_node);
  this->current_time_ = start_node->height();

  model_->resetTime();
  this->searchContemporaries(start_node);
  
  // Skip through model changes
  while ( model_->getNextTime() <= current_time_ ) { 
    model_->increaseTime();
  }

  next();
}


// Sets current_interval_ to the next time interval.
void TimeIntervalIterator::next() {
  if (this->inside_node_ != NULL) {
    this->current_interval_.start_height_ = inside_node_->height();
    this->inside_node_ = NULL;
    return;
  }

  if (current_time_ == DBL_MAX) {
    good_ = false;
    return;
  }

  double start_height = this->current_time_;

  // Ensure that both iterators point into the future to determine the end of the
  // interval 
  if ( start_height >= forest_->model().getNextTime() ) { 
    model_->increaseTime();
  }

  if ( start_height >= node_iterator_.height() ) {
    // Update contemporaries 
    contemporaries()->replaceChildren(*node_iterator_);

    // Pruning
    while ( !(*node_iterator_)->is_last() ) {
      // Prunes the next node BEFORE node_iterator_ gets there, 
      // and does there not invalidate it.
      if (!forest_->pruneNodeIfNeeded((*node_iterator_)->next())) break;
    }

    // Move node iterator forwards
    ++node_iterator_;
  }

  double next_model_change_ = forest_->model().getNextTime();

  assert( current_time_ <= next_model_change_ );
  //std::cout << "current_time: " << current_time_ << " ni_height: " << node_iterator_.height() << std::endl;
  assert( current_time_ <= node_iterator_.height() );

  // Now determine the end of the interval
  if ( node_iterator_.height() <= next_model_change_ ) {
    current_time_ = node_iterator_.height();
  } else  {
    current_time_ = next_model_change_;
  }
  //std::cout << " Next Node: " << node_iterator_.height()
  //          << " Next MC: " << next_model_change_ 
  //          << " CT " << current_time_ << std::endl;

  //Don't return TimeIntervals of length zero, as nothing can happen there...
  if (start_height == current_time_) return next();
  
  this->current_interval_ = TimeInterval(this, 
                                         start_height, 
                                         current_time_);
}

void TimeIntervalIterator::searchContemporariesBottomUp(Node* node, const bool use_buffer) {
  contemporaries()->clear();
  Node* start_node = NULL;

  if ( use_buffer ) {
    assert( node->height() >= contemporaries()->buffer_time() ); 
    // check if the buffered contemporaries are contemporaries of node
    double highest_time = -1;
    for (size_t pop = 0; pop < model()->population_number(); ++pop) {
      auto end = contemporaries()->buffer_end(pop);
      for (auto it = contemporaries()->buffer_begin(pop); it != end; ++it) {
        assert(!(*it)->is_root());
        //std::cout << "Checking " << *it << std::endl;
        // Prune the node if needed
        tmp_child_1_ = (*it);
        tmp_child_2_ = (*it)->first_child();
        while (tmp_child_1_->countChildren() == 1 && forest_->pruneNodeIfNeeded(tmp_child_1_)) {
          tmp_child_1_ = tmp_child_2_;
          if (tmp_child_1_ == NULL ) break;
          tmp_child_2_ = tmp_child_2_->first_child();
        }
        if (tmp_child_1_ == NULL || forest_->pruneNodeIfNeeded(tmp_child_1_)) continue;

        // And add it if it is a contemporary
        if (tmp_child_1_->height() <= node->height() && node->height() < tmp_child_1_->parent_height()) {
          contemporaries()->add(tmp_child_1_);
        }

        // Find the oldest buffered node
        if (tmp_child_1_->height() > highest_time) {
          highest_time = tmp_child_1_->height();
          start_node = tmp_child_1_;
        }
      }
    }
    // The node after the oldest node in the buffer should be the first node
    // above the buffers_height.
    assert( start_node != NULL );
    start_node = start_node->next();
    assert( start_node->height() >= contemporaries()->buffer_time() );
  } else {
    start_node = forest()->nodes()->first(); 
  }

  for (NodeIterator ni = forest_->nodes()->iterator(start_node); *ni != node; ++ni) {
    assert(ni.good());

    // Check if *ni is a contemporary of node 
    if ( (*ni)->parent_height() > node->height() ) {
      // Is is; it may however be a node we need to prune
      if ((*ni)->is_first()) tmp_prev_node_ = NULL;
      else tmp_prev_node_ = (*ni)->previous();
      tmp_child_1_ = (*ni)->first_child();

      if (forest_->pruneNodeIfNeeded(*ni)) {
        // Removing the node invalidates the ni
        if (tmp_prev_node_ == NULL) ni = forest_->nodes()->iterator();
        else ni = forest_->nodes()->iterator(tmp_prev_node_);

        // Maybe a child of the node became a contemporary by removing the node
        // This can only happen if the node has only one child.
        if ( tmp_child_1_ != NULL && tmp_child_1_->parent_height() > node->height() ) { 
          this->contemporaries()->add(tmp_child_1_);
        }
      } else {
        // No pruning => Just add to contemporaries
        this->contemporaries()->add(*ni);
      }
    }
  }
}



// Recalculates the borders of the current interval. 
// Used after one or more nodes where created inside the interval due to events
// occurring within.
void TimeIntervalIterator::recalculateInterval() {
  if (!node_iterator_.good()) {
    node_iterator_ = NodeIterator(forest_->nodes()->last());
  } 
  else {
    // Set node iterator back to the node at the current start_height
    while ( (*node_iterator_)->height() > current_interval_.start_height() ) --node_iterator_;
    ++node_iterator_;
  }

  // Then go to the next node
  current_interval_.end_height_ = (*node_iterator_)->height();
  current_time_ = (*node_iterator_)->height();
}

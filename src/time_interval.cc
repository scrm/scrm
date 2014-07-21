/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013, 2014 Paul R. Staab, Sha (Joe) Zhu and Gerton Lunter
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


// Uniformly samples a random node from the current contemporaries.
// Distribution checked.
Node* TimeInterval::getRandomContemporary(const size_t population) {
  assert( population < forest().model().population_number() );
  assert( this->numberOfContemporaries(population) > 0 );

  // Sample the position of the Node we return
  size_t sample = forest().random_generator()->sampleInt(numberOfContemporaries(population));
  
  for (auto it = contemporaries(population).begin(); it != contemporaries(population).end(); ++it) {
    assert( *it != NULL );
    if ( sample == 0 ) return (*it);
    --sample;
  }

  throw std::logic_error("Could not find the contemporary node I wanted to sample :(");
  return NULL;
}

/* --------------------------------------------------------------------
 * TimeIntervalIterator
 * -------------------------------------------------------------------*/

TimeIntervalIterator::TimeIntervalIterator() {
  this->forest_ = NULL;
  this->node_iterator_ = NULL;
  this->good_ = true;
  this->inside_node_ = NULL;
}


TimeIntervalIterator::TimeIntervalIterator(Forest* forest, 
                                           Node* start_node) {
  this->forest_ = forest;
  this->good_ = true;
  this->inside_node_ = NULL;
  this->node_iterator_ = forest->nodes()->iterator(start_node);
  this->current_time_ = start_node->height();
  forest->writable_model()->resetTime();

  // Save the contemporaries for each population in an unordered_set
  std::vector< std::unordered_set<Node*> >(forest->model().population_number());
  
  this->contemporaries_ = std::vector<std::unordered_set<Node*> >(forest->model().population_number());
  this->searchContemporaries(start_node);
  
  // Skip through model changes
  while ( forest_->model().getNextTime() <= start_node->height() ) { 
    forest_->writable_model()->increaseTime();
  }

  next();
}


// Sets current_interval_ to the next time interval.
void TimeIntervalIterator::next() {
  if (this->inside_node_ != NULL ) {
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
    forest_->writable_model()->increaseTime();
  }

  if ( start_height >= node_iterator_.height() ) {
    updateContemporaries(*node_iterator_);

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


void TimeIntervalIterator::updateContemporaries(Node* current_node) {
  tmp_child_1_ = current_node->first_child();
  tmp_child_2_ = current_node->second_child();

  // Don't add root nodes
  assert( current_node != NULL );
  if (!current_node->is_root()) this->addToContemporaries(current_node);
  if (tmp_child_1_ != NULL) { 
    this->removeFromContemporaries(tmp_child_1_);
    if (tmp_child_2_ != NULL) this->removeFromContemporaries(tmp_child_2_);
  }
}

/** 
 * Finds all nodes which code for branches at the height of a given node in the
 * tree (i.e. the node's contemporaries). Saves this nodes in the contemporaries_
 * member.
 * 
 * @param node Node for which we are searching contemporaries
 * @return Nothing. Nodes are saved in contemporaries_.
 */
void TimeIntervalIterator::searchContemporaries(Node *node) {
  // Prefer top-down search if we have many samples and the node is old.
  // This values are empirically optimized.
  if (forest_->nodes()->size() < 2000 ||
      node->height()< 0.0025 * forest_->model().default_pop_size) {
    searchContemporariesBottomUp(node);
  } else {
    searchContemporariesTopDown(node);
  } 
}

void TimeIntervalIterator::searchContemporariesBottomUp(Node* node) {
  clearContemporaries();

  for (NodeIterator ni = forest_->nodes()->iterator(); *ni != node; ++ni) {
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
          this->addToContemporaries(tmp_child_1_);
        }
      } else {
        // No pruning => Just add to contemporaries
        this->addToContemporaries(*ni);
      }
    }
  }
}

/** 
 * Does the same as searchContemporariesOfNode, but works its way recursively
 * down the tree.
 *
 * @param node Node for which we are searching contemporaries
 * @return Nothing. Nodes are saved in contemporaries_.
 */
void TimeIntervalIterator::searchContemporariesTopDown(Node* node, 
                                                       Node* current_node) {
  
  if (current_node == NULL) {
    assert(forest_->primary_root() != NULL);
    clearContemporaries();
    searchContemporariesTopDown(node, forest_->primary_root());
    for (auto it = forest_->secondary_roots_.begin();
         it != forest_->secondary_roots_.end(); ++it) {
      searchContemporariesTopDown(node, *it);
    }
    return;
  }

  if (node == current_node) return;

  // If the current_node is above node's height => no contemporary
  // This is the most common case
  if (current_node->height() > node->height()) {
    // It is important that second child comes first here, because it might
    // prune itself away.
    if (current_node->second_child() != NULL) searchContemporariesTopDown(node, current_node->second_child());
    if (current_node->first_child() != NULL) searchContemporariesTopDown(node, current_node->first_child());
    return;
  } 

  // If we are here, current_node is a contemporary, unless it needs to be pruned 
  tmp_child_1_ = current_node->first_child();
  tmp_child_2_ = current_node->second_child();

  if (forest_->pruneNodeIfNeeded(current_node, false)) {
    // We just pruned a node that would have been a contemporary.
    // If it had only one child, that might be a contemporary now.
    if (tmp_child_2_ == NULL && tmp_child_1_ != NULL) {
      searchContemporariesTopDown(node, tmp_child_1_);
    } 
    // If the node was old or orphaned, there is nothing we have to do. 
    return;
  }

  // If it was not pruned, it really is a contemporary
  if (!current_node->is_root()) {
    this->addToContemporaries(current_node);
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

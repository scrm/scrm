#include "time_interval.h"

/* --------------------------------------------------------------------
 * TimeInterval
 * -------------------------------------------------------------------*/

TimeInterval::TimeInterval() {
  this->forest_ = NULL;
  this->start_height_ = 0;
  this->end_height_ = 0;
}


TimeInterval::TimeInterval(Forest* forest, 
             double start_height, 
             double end_height,
             std::set<Node*> contemporaries) {

  this->forest_ = forest;
  this->start_height_ = start_height;
  this->end_height_ = end_height;
  this->contemporaries_ = contemporaries;
}


// Uniformly samples a random node from the current contemporaries.
// Distribution checked.
Node* TimeInterval::getRandomContemporary() const {
  if ( this->contemporaries().size() == 0) 
    throw std::out_of_range("Error: Sampling from empty contemporaries");
  size_t sample = this->forest_->random_generator()->sampleInt(this->contemporaries().size());

  std::set<Node*>::iterator it = contemporaries_.begin();
  while (sample > 0) {
    ++it;
    --sample;
  }
  
  return *it;
}


//Removes a Node from the contemporaries if it is there.
void TimeInterval::removeFromContemporaries(Node* node) {
  contemporaries_.erase(node);
}




/* --------------------------------------------------------------------
 * TimeIntervalIterator
 * -------------------------------------------------------------------*/

TimeIntervalIterator::TimeIntervalIterator() {
  this->forest_ = NULL;
  this->node_iterator_ = forest_->nodes()->iterator();
  this->good_ = true;
  this->inside_node_ = NULL;
};


TimeIntervalIterator::TimeIntervalIterator(Forest* forest, Node const* start_node) {
  this->forest_ = forest;
  this->good_ = true;
  this->inside_node_ = NULL;
  this->node_iterator_ = forest->nodes()->iterator();
  
  //Skipt intervals below start_height
  for( ; *node_iterator_ != start_node; node_iterator_++ ) {
    if ( ! node_iterator_.good() )
      throw std::out_of_range("TimeIntervalIterator: start_node not found");

    if ( (*node_iterator_)->parent_height() > start_node->height() )
      contemporaries_.insert(*node_iterator_);
  }
  
  next();
}


void TimeIntervalIterator::next() {
  dout << "1" << std::endl;
  if (this->inside_node_ != NULL ) {
    this->current_event_.start_height_ = inside_node_->height();
    // end_height doesn't change
    // contemporaries are already changed
    this->inside_node_ = NULL;
    return;
  }
  dout << "1.5" << std::endl;

  if (!node_iterator_.good()) {
    good_ = false;
    return;
  }

  dout << "1.75" << std::endl;
  dout << *node_iterator_ << std::endl;
  double start_height = (*node_iterator_)->height();
  dout << "2" << std::endl;

  //Update contemporaries
  if ( (*node_iterator_)->lower_child() != NULL ) 
    contemporaries_.erase((*node_iterator_)->lower_child());
  if ( (*node_iterator_)->higher_child() != NULL ) 
    contemporaries_.erase((*node_iterator_)->higher_child());
  
  if ( !(*node_iterator_)->is_root() ) 
    contemporaries_.insert(*node_iterator_);

  ++node_iterator_;
  double end_height;
  if (node_iterator_.good()) end_height = (*node_iterator_)->height();
  else end_height = FLT_MAX;
 

  //Don't return TimeIntervals of length zero, as nothing can happen there...
  if (start_height == end_height) return next();
  
  this->current_event_ = TimeInterval(this->forest_, start_height, end_height, contemporaries_); 
}

bool TimeIntervalIterator::good() const {
 return this->good_;
}


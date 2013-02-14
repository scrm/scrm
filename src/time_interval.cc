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
             std::set<Node*>* contemporaries) {

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

  std::set<Node*>::iterator it = contemporaries().begin();
  while (sample > 0) {
    ++it;
    --sample;
  }
  
  return *it;
}


bool TimeInterval::checkContemporaries() const {
  for (std::set<Node*>::iterator it = contemporaries().begin(); it != contemporaries().end(); ++it) {
    if ( (*it)->height() > start_height_ || (*it)->parent_height() < end_height_ ) {
      std::cout << "Non-comtemporary node " << *it << " in contemporaries" << std::endl; 
      return 0;
    }
  }
  return 1;
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
  if (this->inside_node_ != NULL ) {
    this->current_event_.start_height_ = inside_node_->height();
    this->inside_node_ = NULL;
    assert( this->current_event_.checkContemporaries() );  
    return;
  }

  if (!node_iterator_.good()) {
    good_ = false;
    return;
  }

  double start_height = (*node_iterator_)->height();

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
  
  this->current_event_ = TimeInterval(this->forest_, start_height, end_height, &contemporaries_);
  assert( this->current_event_.checkContemporaries() );  
}


bool TimeIntervalIterator::good() const {
 return this->good_;
}


//Removes a Node from the contemporaries if it is there.
void TimeIntervalIterator::removeFromContemporaries(Node* node) {
  this->contemporaries_.erase(node);
}


// Recalculates the borders of the current interval. 
// Used after one or more nodes where created inside the interval due to events
// occurring within.
void TimeIntervalIterator::recalculateInterval() {
  // Set node iterator back to the node at the current start_height
  while ( (*node_iterator_)->height() > current_event_.start_height() ) --node_iterator_;
  // Than go to the next node
  ++node_iterator_;
  current_event_.end_height_ = (*node_iterator_)->height();
}

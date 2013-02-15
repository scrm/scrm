#include "time_interval.h"

/* --------------------------------------------------------------------
 * TimeInterval
 * -------------------------------------------------------------------*/

TimeInterval::TimeInterval() {
  this->forest_ = NULL;
  this->start_height_ = 0;
  this->end_height_ = 0;
  this->contemporaries_ = NULL;
}


TimeInterval::TimeInterval(Forest* forest, 
             double start_height, 
             double end_height,
             std::vector<Node*> *contemporaries) {

  assert( contemporaries != NULL );
  this->forest_ = forest;
  this->start_height_ = start_height;
  this->end_height_ = end_height;
  this->contemporaries_ = contemporaries;
}


// Uniformly samples a random node from the current contemporaries.
// Distribution checked.
Node* TimeInterval::getRandomContemporary() const {
  if ( this->numberOfContemporaries() == 0) 
    throw std::out_of_range("Error: Sampling from empty contemporaries");

  size_t sample = this->forest_->random_generator()->sampleInt(this->numberOfContemporaries());
  return contemporaries().at(sample);
}


// Checks if all nodes in contemporaries are contempoaries.
// Does not check if their are an nodes missing.
bool TimeInterval::checkContemporaries() const {
  //std::cout << "Cont: " << contemporaries_ << std::endl;
  assert( contemporaries_ != NULL );
  for (std::vector<Node*>::const_iterator it = contemporaries_->begin(); 
       it != contemporaries_->end(); ++it) {

    //std::cout << "Size " << numberOfContemporaries() << " checking " << *it << std::endl;
    if ( (*it)->height() > start_height_ || (*it)->parent_height() < end_height_ ) {
      std::cout << "Non-contemporary node " << *it << " in contemporaries" << std::endl; 
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
  //std::cout << "New Iterator" << std::endl;
  this->forest_ = forest;
  this->good_ = true;
  this->inside_node_ = NULL;
  this->node_iterator_ = forest->nodes()->iterator();
  this->contemporaries_ = std::vector<Node*>(0);
  
  //Skipt intervals below start_height
  for( ; *node_iterator_ != start_node; node_iterator_++ ) {
    if ( ! node_iterator_.good() )
      throw std::out_of_range("TimeIntervalIterator: start_node not found");

    if ( (*node_iterator_)->parent_height() > start_node->height() )
      this->addToContemporaries(*node_iterator_);
  }
  
  next();
  //std::cout << "New Iterator Created" << std::endl;
}


// Sets current_interval_ to the next time interval.
void TimeIntervalIterator::next() {
  if (this->inside_node_ != NULL ) {
    this->current_interval_.start_height_ = inside_node_->height();
    this->inside_node_ = NULL;
    assert( this->current_interval_.checkContemporaries() );  
    return;
  }

  if (!node_iterator_.good()) {
    good_ = false;
    return;
  }

  double start_height = (*node_iterator_)->height();

  //Update contemporaries
  if ( (*node_iterator_)->lower_child() != NULL ) 
    this->removeFromContemporaries((*node_iterator_)->lower_child());
  if ( (*node_iterator_)->higher_child() != NULL ) 
    this->removeFromContemporaries((*node_iterator_)->higher_child());
  
  if ( !(*node_iterator_)->is_root() ) 
    this->addToContemporaries(*node_iterator_);

  ++node_iterator_;
  double end_height;
  if (node_iterator_.good()) end_height = (*node_iterator_)->height();
  else end_height = FLT_MAX;

  //Don't return TimeIntervals of length zero, as nothing can happen there...
  if (start_height == end_height) return next();
  
  this->current_interval_ = TimeInterval(this->forest_, start_height, end_height, &contemporaries_);
  assert( this->current_interval_.checkContemporaries() );  
}


// Removes a Node from the contemporaries if it is there.
// Currently their are situations where we try to do this also for nodes not in
// contemporaries. Maybe we could optimize this, but this may not be easy.
void TimeIntervalIterator::removeFromContemporaries(Node* node) {
  //std::cout << "Removing " << node << std::endl;
  std::vector<Node*>::iterator it;
  for (it = contemporaries_.begin(); it != contemporaries_.end(); ++it) {
    //std::cout << *it << " == " << node << " ?" << std::endl;
    if ( *it == node ) {
      contemporaries_.erase(it);
      return;
    }
  }
  //throw std::logic_error("Trying to delete noexisting node from contemporaries");
}


// Recalculates the borders of the current interval. 
// Used after one or more nodes where created inside the interval due to events
// occurring within.
void TimeIntervalIterator::recalculateInterval() {
  // Set node iterator back to the node at the current start_height
  while ( (*node_iterator_)->height() > current_interval_.start_height() ) --node_iterator_;
  // Than go to the next node
  ++node_iterator_;
  current_interval_.end_height_ = (*node_iterator_)->height();
}

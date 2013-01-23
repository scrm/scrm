#include "event.h"

Event::Event() {
  this->forest_ = NULL;
  this->start_height_ = 0;
  this->end_height_ = 0;
}

Event::Event(Forest* forest, 
             double start_height, 
             double end_height,
             std::vector<Node*> contemporaries) {

  this->forest_ = forest;
  this->start_height_ = start_height;
  this->end_height_ = end_height;
  this->contemporaries_ = contemporaries;
}

Node* Event::getRandomContemporary() const {
  int sample = this->forest_->random_generator()->sampleInt(this->contemporaries().size());
  return(this->contemporaries()[sample]);
}

//Removes a Node from the contemporaries if it is there.
void Event::removeFromContemporaries(Node* node) {
  for (std::vector<Node*>::iterator it = contemporaries_.begin();
       it!=contemporaries_.end(); ++it) {
    if ( *it == node ) {
      contemporaries_.erase(it);
      return;
    }
  }
  throw std::logic_error("Trying to remove non-existing node");
}



EventIterator::EventIterator() {
  this->forest_ = NULL;
  this->contemporaries_ = std::vector<Node*>(0);
  this->twig_iterator_ = forest_->nodes()->iterator();
  this->good_ = true;
};


EventIterator::EventIterator(Forest* forest, Node const* start_node) {
  this->forest_ = forest;
  this->contemporaries_ = std::vector<Node*>(0);
  this->contemporaries_.reserve(forest_->nodes()->size());
  this->good_ = true;

  //Skipt intervals below start_height
  for( twig_iterator_ = forest->nodes()->iterator();
       *twig_iterator_ != start_node; twig_iterator_++ ) {
    if ( ! twig_iterator_.good() )
      throw std::out_of_range("EventIterator: start_node not found");

    if ( (*twig_iterator_)->parent_height() > start_node->height() )
      contemporaries_.push_back(*twig_iterator_);
  }
  
  next();
}



void EventIterator::next() {
  if (!twig_iterator_.good()) {
    good_ = false;
    return;
  }

  double start_height = (*twig_iterator_)->height();

  if (!( (*twig_iterator_)->is_root() || (*twig_iterator_)->is_fake())) 
    contemporaries_.push_back(*twig_iterator_);

  ++twig_iterator_;
  double end_height;
  if (twig_iterator_.good()) end_height = (*twig_iterator_)->height();
  else end_height = FLT_MAX;
 
  //Don't return Events of length zero, as nothing can happen there...
  if (start_height == end_height) return next();

  //Update contemporaries
  //This could be a bit more optimized
  int i=0;
  while ( i < (contemporaries_.size()) ) {
    if ( contemporaries_[i]->parent_height() <= start_height ) {
      contemporaries_.erase(contemporaries_.begin() + i);
    } else {
      ++i;
    }
  }
  
  this->current_event_ = Event(this->forest_, start_height, end_height, contemporaries_); 
  // TODO: Remove; compatibility only atm
}

bool EventIterator::good() const {
 return this->good_;
}


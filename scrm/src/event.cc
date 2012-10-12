#include "event.h"

Event::Event() {
  this->start_height_ = 0;
  this->end_height_ = 0;
}

Event::Event(double start_height, double end_height, std::vector<Node*> contemporaries){
  this->start_height_ = start_height;
  this->end_height_ = end_height;
  this->contemporaries_ = contemporaries;
}


EventIterator::EventIterator() {
  this->forest_ = NULL;
  this->start_height_ = 0;
  this->contemporaries_ = std::vector<Node*>(0); 
  this->twig_iterator_ = std::vector<Node*>::iterator();
};
  
EventIterator::EventIterator(Forest* forest, const double &start_height) {
  if (start_height < 0) 
    throw std::out_of_range("EventIterator: start_height negative!");
  
  this->forest_ = forest;
  this->twig_iterator_ = forest->getNodeFwIterator();
  this->start_height_ = start_height;
  this->contemporaries_ = std::vector<Node*>(0);
  this->contemporaries_.reserve(forest_->countNodes());

  //Skipt intervals below start_height
  while ((*twig_iterator_)->height() < start_height) {
    if ( twig_iterator_ == forest_->getNodesEnd()) break;
    //Does the twig end above start_height?
    if ((*twig_iterator_)->parent_height() > start_height){
      contemporaries_.push_back(*twig_iterator_);
    }

    ++twig_iterator_;
  }
}

EventIterator::~EventIterator() { };
  
Event EventIterator::next() {
  if (twig_iterator_ == forest_->getNodesEnd())
    throw std::out_of_range("EventIterator out of range");
  
  double start_height = (*twig_iterator_)->height();

  //Return the rest of the interval start_height_ is in (if any)
  if (start_height_ > 0) {
    double cur_height = start_height_;
    start_height_ = 0;
    if ( start_height > cur_height ) {
      return Event(cur_height, start_height, contemporaries_);
    }
  }
  
  contemporaries_.push_back(*twig_iterator_); 
  ++twig_iterator_;

  //Set interval borders
  double end_height;
  if (twig_iterator_ == forest_->getNodesEnd()) 
    end_height = FLT_MAX;
  else
    end_height = (*twig_iterator_)->height();

  //Update contemporaries
  //This could be a bit more optimized
  for (int i=0; i < (contemporaries_.size() - 1); i++) {
    if ( contemporaries_[i]->parent_height() <= start_height ) {
      contemporaries_.erase(contemporaries_.begin() + i);
      --i; 
    }
  }

  //Don't return Events of length zero, as nothing can happen there...
  if (start_height == end_height) return next();
  
  return Event(start_height, end_height, contemporaries_);
}

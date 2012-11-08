#include "event.h"

Event::Event() {
  this->forest_ = NULL;
  this->start_height_ = 0;
  this->end_height_ = 0;
}

Event::Event(Forest* forest, 
             double start_height, 
             double end_height, 
             std::vector<Node*> contemporaries){
  this->forest_ = forest;
  this->start_height_ = start_height;
  this->end_height_ = end_height;
  this->contemporaries_ = contemporaries;
}

Node* Event::getRandomContemporary() {
  int sample = this->forest_->random_generator()->sampleInt(this->contemporaries().size());
  return(this->contemporaries()[sample]);
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
    if ((*twig_iterator_)->height_above() > start_height){
      if (!(*twig_iterator_)->is_root()) 
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
      return Event(this->forest_, cur_height, start_height, contemporaries_);
    }
  }
  
  if (!(*twig_iterator_)->is_root()) contemporaries_.push_back(*twig_iterator_);
  ++twig_iterator_;

  //Set interval borders
  double end_height;
  if (twig_iterator_ == forest_->getNodesEnd()) 
    end_height = FLT_MAX;
  else
    end_height = (*twig_iterator_)->height();
  
  //Don't return Events of length zero, as nothing can happen there...
  if (start_height == end_height) return next();

  //Update contemporaries
  //This could be a bit more optimized
  for (int i=0; i < (contemporaries_.size() - 1); i++) {
    if ( contemporaries_[i]->parent_height() <= start_height ) {
      contemporaries_.erase(contemporaries_.begin() + i);
      --i; 
    }
  }
  
  return Event(this->forest_, start_height, end_height, contemporaries_);
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

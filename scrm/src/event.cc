#include "event.h"

EventIterator::EventIterator() {
  this->forest_ = NULL;
  this->current_height_ = -1;
  this->contemporaries_ = std::vector<Node*>(0); 
  this->twig_iterator_ = std::vector<Node*>::iterator();
};
  
EventIterator::EventIterator(Forest* forest, const double &start_height) {
  this->forest_ = forest;
  this->twig_iterator_ = forest->getNodeFwIterator();
  this->current_height_ = start_height;
  this->contemporaries_ = std::vector<Node*>(forest_->countNodes()); 

  //Skipt intervals below start_height
  while ((*twig_iterator_)->height() < start_height) {
    //Does the twig end above start_height?
    if ((*twig_iterator_)->height() + (*twig_iterator_)->height_above() > start_height){
      contemporaries_.push_back(*twig_iterator_);
    } 
    ++twig_iterator_;
  }
}

EventIterator::~EventIterator() { };
  
struct Event EventIterator::next() {
  if (twig_iterator_ == forest_->getNodesEnd())
    throw std::out_of_range("EventIterator out of range");
  
  double start_height = (*twig_iterator_)->height();

  ++twig_iterator_;

  double end_height;
  if (twig_iterator_ == forest_->getNodesEnd()) 
    end_height = FLT_MAX;
  else
    end_height = (*twig_iterator_)->height();

  Event event = {start_height, end_height, contemporaries_};
  return(event);
}

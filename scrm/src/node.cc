#include "node.h"
   
Node::Node() { init(); };
Node::Node(double height) { init(height); };
Node::Node(double height, bool active) { init(height, active); };
Node::~Node() {};

void Node::init(double height, bool active) {
  this->set_height(height);
  this->set_active(active);
  this->set_parent(NULL);
  this->set_higher_child(NULL); 
  this->set_lower_child(NULL);
}
   
double Node::parent_height() { 
  if ( this->parent() == NULL ) return FLT_MAX;
  return this->parent()->height(); 
}

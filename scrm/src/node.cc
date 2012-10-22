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

void Node::change_child(Node* from, Node* to) {
  if ( this->higher_child() == from )
    this->set_higher_child(to);
  else if ( this->lower_child() == from )
    this->set_lower_child(to);
  else throw std::invalid_argument("Can't find child node to replace");

  if (this->higher_child()->height() < this->lower_child()->height()) {
    Node* tmp = this->higher_child();
    this->set_higher_child(this->lower_child());
    this->set_lower_child(tmp);
  }
}

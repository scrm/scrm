#include "node.h"

Node::Node() { init(); };
Node::Node(double height) { init(height, true); };
Node::Node(double height, bool active) { init(height, active); };
Node::~Node() {};

void Node::init(double height, bool active) {
  this->set_height(height);
  this->set_active(active);
  this->set_parent(NULL);
  this->set_higher_child(NULL); 
  this->set_lower_child(NULL);
}

Node* Node::parent() const {
  if (this->parent_ == NULL) 
    throw std::logic_error("Trying to access parent of ultimate root"); 
  return this->parent_; 
}

double Node::parent_height() const {
  if ( this->is_root() ) return FLT_MAX;
  return this->parent()->height(); 
}

void Node::change_child(Node* from, Node* to) {
  if ( this->higher_child() == from )
    this->set_higher_child(to);
  else if ( this->lower_child() == from )
    this->set_lower_child(to);
  else throw std::invalid_argument("Can't find child node to replace");

  if (this->numberOfChildren() < 2) return;
  if (this->higher_child()->height() < this->lower_child()->height()) {
    Node* tmp = this->higher_child();
    this->set_higher_child(this->lower_child());
    this->set_lower_child(tmp);
  }
}

int Node::numberOfChildren() { 
  return( (this->higher_child() != NULL) + (this->lower_child() != NULL) );
}

//TRUE iff we are the root of a real tree
//false for the ultimate_root
bool Node::is_root() const {
  if ( this->is_fake() ) return false;
  //Special case: toot of a tree just cut away for coalescence. These are not
  //added to the fake tree to increase performance.
  if ( this->parent_ == NULL ) return true;
  //or root of a tree:
  return ( this->parent()->is_fake() );
}
   
bool Node::is_fake() const {
  return ( this->height() == FLT_MAX ); 
}

bool Node::is_ultimate_root() const {
  return ( this->is_fake() && this->parent_ == NULL );
}

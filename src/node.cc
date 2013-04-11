#include "node.h"

Node::Node() { init(); };
Node::Node(double height) { init(height, true); };
Node::Node(double height, bool local) { init(height, local); };
Node::Node(double height, bool local, size_t last_update) { 
  init(height, local, last_update); };
Node::Node(double height, bool local, size_t last_update, size_t samples_below)
  { init(height, local, last_update, samples_below); };
Node::Node(double height, bool local, size_t last_update, size_t samples_below, double length_below)
  { init(height, local, last_update, samples_below, length_below); };
Node::~Node() {};

void Node::init(double height, bool local, size_t last_update,
                size_t samples_below, double length_below) {
  this->set_height(height);
  this->set_local(local);
  this->set_parent(NULL);
  this->set_higher_child(NULL);
  this->set_lower_child(NULL);
  this->set_last_update(last_update);
  this->set_samples_below(samples_below);
  this->set_length_below(length_below);
  this->set_previous(NULL);
  this->set_next(NULL);
}

Node* Node::parent() const {
  //if (this->parent_ == NULL) 
  //  throw std::logic_error("Trying to access parent of ultimate root"); 
  return this->parent_; 
}

void Node::change_child(Node* from, Node* to) {
  if ( this->lower_child() == from )
    this->set_lower_child(to);
  else if ( this->higher_child() == from )
    this->set_higher_child(to);
  else throw std::invalid_argument("Can't find child node to replace");

  if (this->numberOfChildren() < 2) return;
  this->sort_children();
}

void Node::sort_children() {
  if (this->higher_child()->height() < this->lower_child()->height()) {
    Node* tmp = this->higher_child();
    this->set_higher_child(this->lower_child());
    this->set_lower_child(tmp);
  }
}

int Node::numberOfChildren() const { 
  return( (this->higher_child() != NULL) + (this->lower_child() != NULL) );
}

void Node::remove_child(Node* child) {
  if ( this->lower_child() == child ) {
    this->set_lower_child(this->higher_child());
    this->set_higher_child(NULL);
    return;
  } 

  if ( this->higher_child() == child ) {
    this->set_higher_child(NULL);
    return;
  }

  throw std::invalid_argument("Can't find child to delete");
}

bool Node::in_sample() const {
  if ( this->height() == 0 ) return true;
  return false;
}



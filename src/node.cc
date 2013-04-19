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
  this->set_second_child(NULL);
  this->set_first_child(NULL);
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
  if ( this->first_child() == from )
    this->set_first_child(to);
  else if ( this->second_child() == from )
    this->set_second_child(to);
  else throw std::invalid_argument("Can't find child node to replace");
}

int Node::numberOfChildren() const { 
  if (first_child() == NULL) return 0;
  else if (second_child() == NULL) return 1;
  else return 2;
}

void Node::remove_child(Node* child) {
  if ( this->first_child() == child ) {
    this->set_first_child(this->second_child());
    this->set_second_child(NULL);
    return;
  } 

  if ( this->second_child() == child ) {
    this->set_second_child(NULL);
    return;
  }

  throw std::invalid_argument("Can't find child to delete");
}

bool Node::in_sample() const {
  if ( this->height() == 0 ) return true;
  return false;
}



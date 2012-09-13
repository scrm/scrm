#include "node.h"
   
Node::Node() { init(); };
Node::Node(int height) { init(height); };
Node::Node(int height, bool active) { init(height, active); };
Node::~Node() {};

void Node::init(int height, bool active, 
             Node *parent, Node *higher_child, Node *lower_child) {

  this->set_height(height);
  this->set_active(active);
  this->set_parent(parent);
  this->set_higher_child(higher_child);
  this->set_lower_child(lower_child);

}

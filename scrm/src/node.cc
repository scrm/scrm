#include "node.h"
#include <cstddef>
using namespace std;

Node::Node() {
  this->set_height(0);
  this->set_parent(NULL);
  this->set_lower_child(NULL);
}

Node::Node(int height) {
  this->set_height(height);
  this->set_parent(NULL);
  this->set_lower_child(NULL);
}

Node::~Node() {

}

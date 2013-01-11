#include "node_container.h"

/*******************************************************
 * Node Container
 *******************************************************/

NodeContainer::NodeContainer() {
  nodes_.reserve(100);
  unsorted_node_ = NULL;
};

NodeContainer::~NodeContainer() {
  std::vector<Node*>::iterator it;
  for (it = nodes_.begin(); it!=nodes_.end(); ++it) {
    delete *it;
  }

  if (unsorted_node_ != NULL) delete unsorted_node_;
}

NodeIterator NodeContainer::interator() { 
  return( NodeIterator(*this) ); 
}

Node* NodeContainer::get(size_t nr, const bool &from_sorted) {
  //if (from_sorted) 
  return nodes_.at(nr);
  //else return unsorted_nodes_.at(nr);
}

void NodeContainer::add(Node *node, const bool &sort) {
  /*if (!sort) {
    unsorted_nodes_.push_back(node);
    return;
    }*/

  std::vector<Node*>::iterator it;
  for (it = nodes_.begin(); it!=nodes_.end(); ++it) {
    if ((*it)->height() >= node->height()) break;
  }
  nodes_.insert(it, node);
};

void NodeContainer::remove(Node *node) {;
  std::vector<Node*>::iterator it;
  for (it = nodes_.begin(); it!=nodes_.end(); ++it) {
    if ((*it) == node ) break;
  }
  if (it == nodes_.end()) throw std::logic_error("Trying to delete apparently non-existing node");
  nodes_.erase(it);
}

void NodeContainer::move(Node *node, const double &new_height) {};

void NodeContainer::sort() {};
void NodeContainer::is_sorted() {};

size_t NodeContainer::size() {
  return (nodes_.size());
};


/*******************************************************
 * Node Iterator
 *******************************************************/

NodeIterator::NodeIterator() {
  //iter_ = new ;
  //good_ = false;
}

NodeIterator::NodeIterator(const NodeContainer& nc) {
  iter_ = nc.nodes_.begin();
  //good_ = true;
}

bool NodeIterator::good() {
  return( ! ((*iter_)->is_fake()) );
}

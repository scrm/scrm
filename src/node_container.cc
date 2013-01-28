#include "node_container.h"

/*******************************************************
 * Con- & Destructor
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



/*******************************************************
 * Iterators
 *******************************************************/

NodeIterator NodeContainer::iterator() { 
  return(NodeIterator(*this)); 
}

ConstNodeIterator NodeContainer::iterator() const { 
  return(ConstNodeIterator(*this)); 
}

ReverseConstNodeIterator NodeContainer::reverse_iterator() const { 
  return(ReverseConstNodeIterator(*this)); 
}



/*******************************************************
 * Management of Nodes
 *******************************************************/

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
  assert( this->sorted() );
};

void NodeContainer::remove(Node *node) {;
  std::vector<Node*>::iterator it;
  for (it = nodes_.begin(); it!=nodes_.end(); ++it) {
    if ((*it) == node ) break;
  }
  if (it == nodes_.end()) throw std::logic_error("Trying to delete apparently non-existing node");
  nodes_.erase(it);
  assert( this->sorted() );
}

void NodeContainer::move(Node *node, const double &new_height) {
  //Optimize!
  this->remove(node);
  node->set_height(new_height);
  this->add(node);
  assert( this->sorted() );
};

size_t NodeContainer::size() const {
  return nodes_.size();
};



/*******************************************************
 * Consistency checking
 *******************************************************/

bool NodeContainer::sorted() const { 
  double cur_height = 0;
  for (ConstNodeIterator it = iterator(); it.good(); ++it) {
    if ((*it)->height() < cur_height) {
      return(0);
    } else {
      cur_height = (*it)->height();
    }
  }
  return(1);
};

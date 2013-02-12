#ifndef scrm_src_node_container
#define scrm_src_node_container

#include <vector>
#include <stdexcept>
#include <cfloat>
#include <cassert>

#include "node.h"

class NodeIterator;
class ConstNodeIterator;
class ReverseConstNodeIterator;

class NodeContainer {
 public:
  NodeContainer();
  ~NodeContainer();
  
  NodeIterator iterator();
  ConstNodeIterator iterator() const;
  ReverseConstNodeIterator reverse_iterator() const;

  void add(Node *node, const bool &sort=true);
  void remove(Node *node);
  void move(Node *node, const double &new_height);
  void clear() { };

  Node const* get(size_t nr) const { return get(nr); };
  Node* get(size_t nr, const bool &from_sorted=true);
  
  size_t size() const;  
  bool sorted() const; 

#ifdef UNITTEST
  friend class TestNodeContainer;
#endif
  friend class NodeIterator;
  friend class ConstNodeIterator;
  friend class ReverseConstNodeIterator;

 private:
  //const std::vector<Node*> nodes() const { return nodes_; }

  void add_after(Node* add, Node* after);

  Node* first_node_;
  Node* last_node_;

  Node* first() const { return first_node_; };
  Node* last() const { return last_node_; };
  
  void set_first(Node* node) { first_node_ = node; }
  void set_last(Node* node) { last_node_ = node; }

  Node* unsorted_node_;
};


class NodeIterator {
 public:
  NodeIterator() { current_node_ = NULL; };
  NodeIterator(NodeContainer& nc) { current_node_ = nc.first(); };
  ~NodeIterator() {};

  Node* operator*() { return current_node_; };
  
  Node const* operator++() {
    if (current_node_ == NULL) 
      throw std::out_of_range("Node iterator out of range");

    current_node_ = current_node_->next();
    return current_node_;
  }
  
  Node const* operator++(int) {
    if (current_node_ == NULL) 
      throw std::out_of_range("Node iterator out of range");
    
    Node* ret = current_node_;
    current_node_ = current_node_->next();
    return ret;
  }

  bool good() const { 
    return current_node_ != NULL;
  }
  
#ifdef UNITTEST
  friend class TestNodeContainer;
#endif

 private:
  Node* current_node_;
  NodeContainer* nc_;
};


class ConstNodeIterator {
 public:
  ConstNodeIterator() {};
  ConstNodeIterator(const NodeContainer& nc) { 
  }
  ~ConstNodeIterator() {};

  Node const* operator*() { return NULL; }
  
  Node const* operator++() { return NULL;
  }
  
  Node const* operator++(int) { 
    return NULL;
  }
 
  bool good() const { return 1; }

 private:
};


class ReverseConstNodeIterator {
 public:
  ReverseConstNodeIterator() {};
  ReverseConstNodeIterator(const NodeContainer &nc) { 
  };
  ~ReverseConstNodeIterator() {};

  Node const* operator*() { return NULL; }
  Node const* operator++() { return NULL;
  }
  
  Node const* operator++(int) { return NULL;
  }

  bool good() const { return 1;  }

 private:
}; 
#endif

#ifndef scrm_src_node_container
#define scrm_src_node_container

#include <vector>
#include <stdexcept>
#include <cfloat>
#include <cassert>
#include <iostream>

#include "node.h"
#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif
class NodeIterator;
class ConstNodeIterator;
class ReverseConstNodeIterator;

class NodeContainer {
 public:
  NodeContainer();
  //~NodeContainer(){ clear(); }
    ~NodeContainer();
    
  NodeIterator iterator();
  ConstNodeIterator iterator() const;
  ReverseConstNodeIterator reverse_iterator() const;

  void add(Node* node, Node* after_node=NULL);
  void remove(Node *node, const bool &del=true);
  void move(Node *node, const double &new_height);
  void clear();

  Node* at(size_t nr) const; 
  Node const* get(size_t nr) const { return at(nr); }; 
  Node * get_copy(size_t nr)  { return at(nr); }; 

  
  Node* first() const { return first_node_; };
  Node* last() const { return last_node_; };
  
  size_t size() const;  
  bool sorted() const; 
  void print() const;

#ifdef UNITTEST
  friend class TestNodeContainer;
#endif
  friend class NodeIterator;
  friend class ConstNodeIterator;
  friend class ReverseConstNodeIterator;

 private:
  //const std::vector<Node*> nodes() const { return nodes_; }

  void add_before(Node* add, Node* next_node);

  Node* first_node_;
  Node* last_node_;
  
  void set_first(Node* node) { first_node_ = node; }
  void set_last(Node* node) { last_node_ = node; }

  Node* unsorted_node_;
};


class NodeIterator {
 public:
  NodeIterator() { current_node_ = NULL; };
  NodeIterator(NodeContainer& nc) { current_node_ = nc.first(); };
  ~NodeIterator() {};

  Node* operator*() { 
    if (current_node_ == NULL) 
      throw std::out_of_range("Node iterator out of range");
    return current_node_; 
  }
  
  Node* operator++() {
    if (current_node_ == NULL) 
      throw std::out_of_range("Node iterator out of range");

    if ( current_node_->is_last() ) current_node_ = NULL;
    else current_node_ = current_node_->next();
    return current_node_;
  }

  Node* operator--() {
    if (current_node_ == NULL) 
      throw std::out_of_range("Node iterator out of range");

    if ( current_node_->is_first() ) current_node_ = NULL;
    else current_node_ = current_node_->previous();
    return current_node_;
  }
  
  Node* operator++(int) {
    if (current_node_ == NULL) 
      throw std::out_of_range("Node iterator out of range");
    
    Node* ret = current_node_;
    if ( current_node_->is_last() ) current_node_ = NULL;
    else current_node_ = current_node_->next();
    return ret;
  }

  bool good() const { 
    return current_node_ != NULL;
  }

  double height() const {
    if ( good() ) return current_node_->height();
    else return FLT_MAX;
  }
  
#ifdef UNITTEST
  friend class TestNodeContainer;
#endif

 private:
  Node* current_node_;
};


class ConstNodeIterator {
 public:
  ConstNodeIterator() { current_node_ = NULL; };
  ConstNodeIterator(const NodeContainer& nc) { current_node_ = nc.first(); };
  ~ConstNodeIterator() {};

  Node const* operator*() { 
    if (current_node_ == NULL) 
      throw std::out_of_range("Node iterator out of range");
    return current_node_; 
  }
  
  Node const* operator++() { 
    if (current_node_ == NULL) 
      throw std::out_of_range("Node iterator out of range");

    if ( current_node_->is_last() ) current_node_ = NULL;
    else current_node_ = current_node_->next();

    return current_node_;
  }
     
  Node const* operator++(int) { 
    if (current_node_ == NULL) 
      throw std::out_of_range("Node iterator out of range");
    
    Node* ret = current_node_;
    if ( current_node_->is_last() ) current_node_ = NULL;
    else current_node_ = current_node_->next();
    return ret;
  }
 
  bool good() const { 
    return current_node_ != NULL;
  }

 private:
  Node* current_node_;
};


class ReverseConstNodeIterator {
 public:
  ReverseConstNodeIterator() { current_node_ = NULL; };
  ReverseConstNodeIterator(const NodeContainer &nc) {  current_node_ = nc.last(); };
  ~ReverseConstNodeIterator() {};

  Node const* operator*() { 
    if (current_node_ == NULL) 
      throw std::out_of_range("Node iterator out of range");
    return current_node_; 
  }

  Node const* operator++() { 
    if (current_node_ == NULL) 
      throw std::out_of_range("Node iterator out of range");

    if ( current_node_->is_first() ) current_node_ = NULL;
    else current_node_ = current_node_->previous();
    return current_node_;
  }

  Node const* operator++(int) { 
    if (current_node_ == NULL) 
      throw std::out_of_range("Node iterator out of range");
    
    Node* ret = current_node_;
    if ( current_node_->is_first() ) current_node_ = NULL;
    else current_node_ = current_node_->previous();
    return ret;
  }

  bool good() const { 
    return current_node_ != NULL;
  }

 private:
  Node* current_node_;
}; 
#endif

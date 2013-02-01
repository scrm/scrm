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
  void clear() { nodes_.clear(); };

  Node const* get(size_t nr) const { return nodes_.at(nr); }
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
  const std::vector<Node*> nodes() const { return nodes_; }

  std::vector<Node*> nodes_;
  Node* unsorted_node_;
};


class NodeIterator {
 public:
  NodeIterator() {};
  NodeIterator(NodeContainer& nc) { 
    iter_ = nc.nodes_.begin();
    nc_ = &nc;
  }
  ~NodeIterator() {};

  Node* operator*() { return *iter_; }
  Node const* operator++() {
    if ( !good() ) throw std::out_of_range("NodeIterator out of range");
    return *(++iter_); 
  }
  
  Node const* operator++(int) { 
    if ( !good() ) throw std::out_of_range("NodeIterator out of range");
    return *(iter_++); 
  }

  bool good() const { return( iter_ != nc_->nodes_.end() ); }
  
#ifdef UNITTEST
  friend class TestNodeContainer;
#endif

 private:
  std::vector<Node*>::iterator iter_;
  NodeContainer* nc_;
};


class ConstNodeIterator {
 public:
  ConstNodeIterator() {};
  ConstNodeIterator(const NodeContainer& nc) { 
    iter_ = nc.nodes_.begin();
    end_  = nc.nodes_.end(); 
  }
  ~ConstNodeIterator() {};

  Node const* operator*() { return *iter_; }
  
  Node const* operator++() {
    if ( !good() ) throw std::out_of_range("NodeIterator out of range");
    return *(++iter_); 
  }
  
  Node const* operator++(int) { 
    if ( !good() ) throw std::out_of_range("NodeIterator out of range");
    return *(iter_++); 
  }
 
  bool good() const { return( iter_ != end_ ); }

 private:
  std::vector<Node*>::const_iterator iter_;
  std::vector<Node*>::const_iterator end_;
};


class ReverseConstNodeIterator {
 public:
  ReverseConstNodeIterator() {};
  ReverseConstNodeIterator(const NodeContainer &nc) { 
    iter_ = nc.nodes_.rbegin(); 
    end_  = nc.nodes_.rend(); 
  };
  ~ReverseConstNodeIterator() {};

  Node const* operator*() { return *iter_; }
  Node const* operator++() {
    if ( !good() ) throw std::out_of_range("NodeIterator out of range");
    return *(++iter_); 
  }
  
  Node const* operator++(int) { 
    if ( !good() ) throw std::out_of_range("NodeIterator out of range");
    return *(iter_++); 
  }

  bool good() const { return( iter_ != end_ ); }

 private:
  std::vector<Node*>::const_reverse_iterator iter_;
  std::vector<Node*>::const_reverse_iterator end_;
}; 
#endif

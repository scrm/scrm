#ifndef scrm_src_node_container
#define scrm_src_node_container

#include "forest.h"
//class Forest;
class NodeIterator;

class NodeContainer {
 public:
  NodeContainer();
  ~NodeContainer();
  
  NodeIterator interator();

  void add(Node *node, const bool &sort=true);
  void remove(Node *node);
  void move(Node *node, const double &new_height);

  size_t size();  

  void sort();
  void is_sorted();

#ifdef UNITTEST
  friend class TestNodeContainer;
#endif
  friend class NodeIterator;

 private:
  Node* get(size_t nr, const bool &from_sorted=true);

  std::vector<Node*> nodes_;
  Node* unsorted_node_;
};

class NodeIterator {
 public:
  NodeIterator();
  NodeIterator(const NodeContainer& nc);
  ~NodeIterator() {};

  Node* operator*() const { return *iter_; }
  Node* operator++() { return *(++iter_); }
  Node* operator++(int) { return *(iter_++); }

  bool good();
  
#ifdef UNITTEST
  friend class TestNodeContainer;
#endif

 private:
  std::vector<Node*>::const_iterator iter_;
//  bool good_;
};

#endif

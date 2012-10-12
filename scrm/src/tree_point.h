#ifndef scrm_src_treepoint
#define scrm_src_treepoint

#include <cassert>
#include "node.h"

class TreePoint {
public:
  TreePoint() {};
  TreePoint(Node* above_node, double height_above);
  ~TreePoint () {};

  Node* above_node() { return above_node_; }
  double height_above() { return height_above_; }

private:
  Node* above_node_;
  double height_above_;
};


#endif

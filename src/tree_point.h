#ifndef scrm_src_treepoint
#define scrm_src_treepoint

#include <cassert>
#include "node.h"

class TreePoint {
public:
  TreePoint() {};
  TreePoint(Node* base_node, double height, bool relative);
  ~TreePoint () {};

  Node*  base_node()       const { return base_node_; }
  double relative_height() const { return relative_height_; }
  double height()          const { return height_; }

private:
  Node*  base_node_;
  double relative_height_;
  double height_;
};

#endif

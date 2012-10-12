#include "tree_point.h"

TreePoint::TreePoint(Node* above_node, double height_above) {
  assert( height_above > 0 );
  assert( height_above <= above_node->height_above() );
  assert( above_node != NULL );

  above_node_ = above_node;
  height_above_ = height_above;
}

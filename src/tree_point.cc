#include "tree_point.h"

TreePoint::TreePoint(Node* base_node, double height, bool relative) {
  assert( base_node != NULL );
  assert( height >= 0 );

  base_node_ = base_node;
  if (relative) {
    relative_height_ = height;
    height_ = base_node->height() + relative_height_;
  } else {
    relative_height_ = height - base_node->height();
    height_ = height;
  }

  assert( relative_height_ <= base_node->height_above() );
  assert( relative_height_ >= 0 ); 
}


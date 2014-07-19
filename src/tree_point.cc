/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013, 2014 Paul R. Staab, Sha (Joe) Zhu, Dirk Metzler and Gerton Lunter
 * 
 * This file is part of scrm.
 * 
 * scrm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

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
  
  //std::cout << base_node << ": " << base_node->height() 
  //                       << " - " << base_node->parent_height() << std::endl;
  assert( relative_height_ <= base_node->height_above() );
  assert( relative_height_ >= 0 ); 
}


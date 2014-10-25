/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013, 2014 Paul R. Staab, Sha (Joe) Zhu and Gerton Lunter
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

#ifndef scrm_src_treepoint
#define scrm_src_treepoint

#include "macros.h" // Needs to be before cassert

#include <cassert>
#include "node.h"

class TreePoint {
public:
  TreePoint() {};
  TreePoint(Node* base_node, double height, bool relative);

  Node*  base_node()       const { return base_node_; }
  double relative_height() const { return relative_height_; }
  double height()          const { return height_; }

private:
  Node*  base_node_;
  double relative_height_;
  double height_;
};

#endif

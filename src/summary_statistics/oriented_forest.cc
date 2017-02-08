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

#include "oriented_forest.h"

void OrientedForest::calculate(const Forest &forest) {
  segment_length_ = forest.calcSegmentLength();
  if (segment_length_ == 0.0) return;
  has_rec_ = forest.model().has_recombination();

  size_t pos = 2*forest.sample_size()-2;
  generateTreeData(forest.local_root(), pos, 0, forest.model().scaling_factor()); 
}


void OrientedForest::printSegmentOutput(std::ostream &output) const {
  if (segment_length_ == 0.0) return;
  output << "{" ;
  if (has_rec_) output << "\"length\":" << segment_length_ << ", ";

  // Print parents
  output << "\"parents\":[" ;
  for (int parent : parents_) { 
    output << parent << ( parent != 0 ? "," : "" );
  }
  output << "], ";

  // Print heights
  output << "\"node_times\":[" ;
  double tmrca = heights_.back();
  for (double height : heights_) {
    output << height << ( height != tmrca ? "," : "" );
  }
  output << "]}" << std::endl;
}


void OrientedForest::generateTreeData(Node const* node, size_t &pos, int parent_pos, const double scaling_factor) {
  // Samples have a fixed position in the arrays, given by their label.
  if (node->in_sample()) {
    heights_.at(node->label()-1) = node->height() * scaling_factor;
    parents_.at(node->label()-1) = parent_pos;
    return;
  }

  // Otherwise take the position given by pos and decrease it.
  heights_.at(pos) = node->height() * scaling_factor;
  parents_.at(pos) = parent_pos;
  parent_pos = pos--;

  Node* local_child_1 = node->getLocalChild1();

  if (local_child_1 != NULL) {
    Node* local_child_2 = node->getLocalChild2();

    if (local_child_2 != NULL) {
      // Ensure that identical topologies lead to identical labels of nodes
      if (local_child_2->height() > local_child_1->height()) {
        Node* tmp = local_child_1;
        local_child_1 = local_child_2;
        local_child_2 = tmp;
      }
      generateTreeData(local_child_2, pos, parent_pos+1, scaling_factor);
    }

    generateTreeData(local_child_1, pos, parent_pos+1, scaling_factor);
  }
}

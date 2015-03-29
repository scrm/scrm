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
  size_t pos = 2*forest.sample_size()-2;
  generateTreeData(forest.local_root(), pos, 0); 

  output_buffer_ << "{" ;

  if (forest.model().recombination_rate() > 0.0) {
    output_buffer_ << "\"length\":" << forest.calcSegmentLength() << ", " ;
  } 

  // Print parents
  output_buffer_ << "\"parents\":[" ;
  for (int parent : parents_) { 
    output_buffer_ << parent << ( parent != 0 ? "," : "" );
  }
  output_buffer_ << "], ";

  // Print heights
  output_buffer_ << "\"node_times\":[" ;
  double tmrca = heights_.at(2*forest.sample_size()-2);
  for (double height : heights_) {
    output_buffer_ << height * forest.model().scaling_factor()
                   << ( height != tmrca ? "," : "" );
  }
  output_buffer_ << "]}" << std::endl;
}

void OrientedForest::printLocusOutput(std::ostream &output) const {
  output << output_buffer_.str();  
}

void OrientedForest::generateTreeData(Node const* node, size_t &pos, int parent_pos) {
  // Samples have a fixed position in the arrays, given by their label.
  if (node->in_sample()) {
    heights_.at(node->label()-1) = node->height();
    parents_.at(node->label()-1) = parent_pos;
    return;
  }

  // Otherwise take the position given by pos and decrease it.
  heights_.at(pos) = node->height();
  parents_.at(pos) = parent_pos;
  parent_pos = pos--;

  if (node->getLocalChild1() != NULL) {
    generateTreeData(node->getLocalChild1(), pos, parent_pos+1);
    if (node->getLocalChild2() != NULL) {
      generateTreeData(node->getLocalChild2(), pos, parent_pos+1);
    }
  }
}

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

#include "newick_tree.h"
NewickTree::NewickTree(const size_t sample_size) {
  labels_.reserve(sample_size);
  for (size_t i=1; i <= sample_size; ++i) {
    labels_.push_back(std::to_string(i));
  } 
};

void NewickTree::calculate(const Forest &forest) {
  if (forest.model().recombination_rate() == 0.0) {
    output_buffer_ << generateTree(forest.local_root(), forest) 
                   << ";" << std::endl;  
  } else {
    if (forest.calcSegmentLength(forest.model().finite_sites()) == 0.0) return;
    output_buffer_ << "[" << forest.calcSegmentLength(forest.model().finite_sites()) << "]" 
                   << generateTree(forest.local_root(), forest)
                   << ";" << std::endl;
  }
}

void NewickTree::printLocusOutput(std::ostream &output) {
  output << output_buffer_.str();  
  output_buffer_.str("");
  output_buffer_.clear();
}

/**
 * @brief Prints a part of the tree in newick format
 *
 * @param node The root of the subtree that will be printed
 *
 * @return A part of the tree in newick format
 */
std::string NewickTree::generateTree(Node *node, const Forest &forest) {
  if (node->in_sample()) return labels_.at(node->label()-1);

  Node *left = forest.trackLocalNode(node->first_child());
  Node *right = forest.trackLocalNode(node->second_child());

  return "(" + generateTree(left, forest) + ":" + 
         std::to_string((node->height() - left->height()) * forest.model().scaling_factor()) +
         "," + generateTree(right, forest) + ":" + 
         std::to_string((node->height() - right->height()) * forest.model().scaling_factor()) + ")";
}

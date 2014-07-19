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

#include "newick_tree.h"


void NewickTree::calculate(const Forest &forest) {
  segment_length_ = forest.calcSegmentLength();
  if (segment_length_ > 0.0) tree_ = generateTree(forest.local_root(), forest, has_rec_); 
}


void NewickTree::printSegmentOutput(std::ostream &output) const {
  if (segment_length_ == 0.0) return;
  if (has_rec_) {
    double intpart; // dummy variable for modf
    if (modf(segment_length_, &intpart) == 0.0) output << "[" << (size_t)segment_length_ << "]";
    else output << "[" << segment_length_ << "]";
  }
  output << tree_ << ";" << std::endl;
}


/**
 * @brief Prints a part of the tree in newick format
 *
 * @param node The root of the subtree that will be printed
 *
 * @return A part of the tree in newick format
 */
std::string NewickTree::generateTree(Node const* node, const Forest &forest, const bool use_buffer) {
  // Use tree from buffer if possible
  std::map<Node const*, NewickBuffer>::iterator it = buffer_.find(node);
  if (use_buffer && it != buffer_.end()) {
    if (it->second.recombination > node->last_change()) {
      // Check that the buffered tree is correct.
      assert(it->second.tree.compare(generateTree(node, forest, false)) == 0);
      return it->second.tree;
    }
  }

  // Generate a new tree
  std::stringstream tree;
  tree.precision(this->precision_);
  tree.exceptions(std::ios::failbit); 

  if (node->in_sample()) tree << node->label();
  else { 
    Node *left = node->getLocalChild1();
    Node *right = node->getLocalChild2();

    tree << "(" << generateTree(left, forest, use_buffer) << ":" <<
           (node->height() - left->height()) * forest.model().scaling_factor() <<
           "," << generateTree(right, forest, use_buffer) << ":" <<
           (node->height() - right->height()) * forest.model().scaling_factor() << ")";

    // And add to to the buffer
    if (use_buffer) {
      NewickBuffer buf = {forest.current_rec(), tree.str()};
      buffer_[node] = buf; 
    }
  }

  return tree.str();
}


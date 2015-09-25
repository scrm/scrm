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

#ifndef scrm_src_summary_statistic_newick_tree
#define scrm_src_summary_statistic_newick_tree

#include <sstream>
#include <iostream>
#include <string>
#include <map>

#include "summary_statistic.h"
#include "../forest.h"


/**
 * @brief Save buffered trees along with the recombination number at which
 * they where created. 
 */
struct NewickBuffer {
  size_t recombination;  ///< The recombination at which the subtree was created.
  std::string tree;      ///< The subtree itself.
};

class NewickTree : public SummaryStatistic
{
 public:
  NewickTree() : NewickTree(6, true) { }
  NewickTree(size_t precision) : NewickTree(precision, true) { }
  NewickTree(size_t precision, bool has_recombination) { 
    precision_ = precision; 
    has_rec_ = has_recombination;
  }

  ~NewickTree() {}

  //Virtual methods
  void calculate(const Forest &forest);
  void printSegmentOutput(std::ostream &output) const;

  NewickTree* clone() const { return new NewickTree(precision_, has_rec_); };

  void clear() {
    buffer_.clear();
  }

 private:
  std::string generateTree(Node const* node, const Forest &forest, const bool use_buffer);
  std::string tree_;
  double segment_length_;
  size_t precision_;
  bool has_rec_;

  /**
   * A map to buffer already created subtrees indexed by their 
   * root.
   */
  std::map<Node const*, NewickBuffer> buffer_;
};

#endif

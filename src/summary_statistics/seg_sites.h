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

#ifndef scrm_src_summary_statistic_seg_sites
#define scrm_src_summary_statistic_seg_sites

#include <sstream>
#include <iostream>
#include <valarray>

#include "summary_statistic.h"
#include "../forest.h"
#include "../tree_point.h"

class SegSites : public SummaryStatistic
{
 public:
  SegSites() { set_position(0.0); }
  ~SegSites() {}

#ifdef UNITTEST
  friend class TestSummaryStatistics;
#endif

  //Virtual methods
  void calculate(const Forest &forest);
  void printLocusOutput(std::ostream &output) const;
  SegSites* clone() const { return new SegSites(*this); }

  void clear() { 
    positions_.clear();
    haplotypes_.clear();  
    set_position(0.0);
  };

  size_t countMutations() const { return positions_.size(); };

  double position() const { return position_; };
  std::vector<double> const* positions() const { return &positions_; };

  std::valarray<bool> const* getHaplotype(const size_t mutation) const {
    return &(haplotypes_.at(mutation));
  }

 private:
  std::valarray<bool> getHaplotypes(TreePoint mutation, const Forest &Forest); 

  std::vector<double> positions_;
  std::vector<std::valarray<bool>> haplotypes_;	
  void traversal(Node const* node, std::valarray<bool> &haplotype) const;

  void set_position(const double position) { position_ = position; };
  double position_;
};

#endif

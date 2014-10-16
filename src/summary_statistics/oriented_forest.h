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

#ifndef scrm_src_summary_statistic_oriented_forest
#define scrm_src_summary_statistic_oriented_forest

#include <sstream>
#include <iostream>
#include <vector>

#include "summary_statistic.h"
#include "../forest.h"

class OrientedForest : public SummaryStatistic
{
 public:
   OrientedForest(const size_t sample_size) {
    parents_ = std::vector<int>(2*sample_size-1, 0);
    heights_ = std::vector<double>(2*sample_size-1, 0.0);
   }
   ~OrientedForest() {}

   //Virtual methods
   void calculate(const Forest &forest);
   void printLocusOutput(std::ostream &output) const;
   void clear() {
     output_buffer_.str("");
     output_buffer_.clear();
   }

   OrientedForest* clone() const {
     return new OrientedForest(this->parents_.size());
   }

#ifdef UNITTEST
   friend class TestSummaryStatistics;
#endif

 private:
   OrientedForest() {}
   void generateTreeData(Node const* node, size_t &pos, int parent_pos);

   std::vector<int> parents_;
   std::vector<double> heights_;
   std::ostringstream output_buffer_;
};

#endif

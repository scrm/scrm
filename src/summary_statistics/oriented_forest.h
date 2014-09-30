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
   OrientedForest() {};
   ~OrientedForest() {};

   //Virtual methods
   void calculate(const Forest &forest);
   void printSegmentOutput(std::ostream &output) { (void)output; }
   void printLocusOutput(std::ostream &output);

 private:
   void init_OF_label( const Forest &forest );
   void update_OF_label( const Forest &forest );
   void generateTree( const Forest &forest );

   std::vector<size_t> OF_labels;
   std::vector<double> heights;
   size_t tmp_label;
   std::ostringstream output_buffer_;
};

#endif

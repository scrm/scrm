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

#ifndef scrm_src_summary_statistic_fst
#define scrm_src_summary_statistic_fst

#include <sstream>
#include <iostream>
#include <vector>

#include "summary_statistic.h"
#include "../forest.h"
#include "seg_sites.h"
#include "../model.h"

class Fst : public SummaryStatistic
{
 public:
  Fst(SegSites* seg_sites, const Model &model) : seg_sites_(seg_sites) {
     std::cout << " ok, initialize fst "<< std::endl;
     sample_size_ = std::vector<size_t>(model.sample_size());
   }
   ~Fst(){};

   Fst(const Fst &fst) : seg_sites_(fst.seg_sites_) { }

   //Virtual methods
   void calculate(const Forest &forest);
   void printLocusOutput(std::ostream &output) const;
   void clear() { }
   Fst* clone() const { return new Fst(*this); }

 private:
   SegSites* const seg_sites_;
   std::vector<size_t> sample_size_;
   double fst_;
};

#endif

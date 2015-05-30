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

#ifndef scrm_src_summary_statistic_frequency_spectrum
#define scrm_src_summary_statistic_frequency_spectrum

#include <sstream>
#include <iostream>
#include <vector>

#include "../macros.h"
#include <cassert>

#include "summary_statistic.h"
#include "seg_sites.h"
#include "../model.h"
#include "../forest.h"

class FrequencySpectrum : public SummaryStatistic
{
 public:
   FrequencySpectrum(std::shared_ptr<SegSites> seg_sites, const Model &model) : seg_sites_(seg_sites) {
     sfs_ = std::vector<size_t>(model.sample_size() - 1, 0);
     at_mutation_ = 0;
     //total_sfs_ = std::vector<size_t>(model_.sample_size() - 1, 0);
   }

   FrequencySpectrum(const FrequencySpectrum &sp) : seg_sites_(sp.seg_sites_) { }

   //Virtual methods
   void calculate(const Forest &forest);
   void printLocusOutput(std::ostream &output) const;
   void clear() { 
     for (size_t i = 0; i < sfs_.size(); ++i) sfs_.at(i) = 0;
     at_mutation_ = 0;
   }
   FrequencySpectrum* clone() const { return new FrequencySpectrum(*this); }
   std::vector<size_t> const & sfs() const { return sfs_; }

 private:
   std::shared_ptr<SegSites> seg_sites_;
   std::vector<size_t> sfs_;
   //std::vector<size_t> total_sfs_;
   size_t at_mutation_;
};

#endif

/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013 Paul R. Staab, Sha (Joe) Zhu and Gerton Lunter
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

#ifndef scrm_seg
#define scrm_seg

#include "param.h"
#include "forest.h"
#include <valarray>

class SegDataBlock{
 public:
  SegDataBlock();
  ~SegDataBlock(){};

  //SegDataBlock(Forest * forest, double max_length, int max_num_mut);
  SegDataBlock(Forest * forest, int max_num_mut);
  vector <double> positions;
  vector < valarray<int> > haplotypes;	
};


class SegDataContainer {
 public:
  SegDataContainer();
  SegDataContainer(Param const* user_para, Forest const* forest);
  ~SegDataContainer();

  friend std::ostream& operator<<(std::ostream& stream, const SegDataContainer &sdc);

  bool seg_bool() const { return user_para_->seg_bool(); };
  int numseg;

  double total_seq_length() const { return forest_->model().loci_length(); };
  int nsam() const { return forest_->model().sample_size(); };

  void append_new_seg_data(Forest *forest);
  vector<int> calcSiteFrequencies() const;

 private:
  vector <SegDataBlock*> seg_datas_;

  Param const* user_para_;
  Forest const* forest_; // do we need this here? 

  size_t remaining_max_num_mut_;
};


#endif

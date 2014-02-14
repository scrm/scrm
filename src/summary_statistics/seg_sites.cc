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

#include "seg_sites.h"

void SegSites::calculate(const Forest &forest) {
  double position_at = forest.current_base();
  position_at += forest.random_generator()->sampleExpo(forest.local_tree_length() * forest.model().mutation_rate());

  while (position_at < forest.next_base()) {
    TreePoint mutation = forest.samplePoint();
    haplotypes_.push_back(getHaplotypes(mutation, forest));
    positions_.push_back(position_at / forest.model().loci_length());
    position_at += forest.random_generator()->sampleExpo(forest.local_tree_length() * forest.model().mutation_rate());
  }	
}

void SegSites::printLocusOutput(std::ostream &output) {
  output << "segsites: "<< countMutations() << std::endl;
  if ( countMutations() == 0 ) return;

  output << "positions: " << positions_ << std::endl;

  for (size_t i = 0; i < haplotypes_.at(0).size(); i++){
    for (size_t j = 0; j < haplotypes_.size(); j++){
      output << haplotypes_[j][i];
    }
    output <<"\n";
  }

  positions_.clear();
  haplotypes_.clear();
};

std::valarray<bool> SegSites::getHaplotypes(TreePoint mutation, const Forest &forest) {
  std::valarray<bool> haplotype(forest.model().sample_size());
  forest.traversal(mutation.base_node(), haplotype);
  return haplotype;
}

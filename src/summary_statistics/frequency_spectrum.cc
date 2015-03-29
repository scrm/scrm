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

#include "frequency_spectrum.h"

void FrequencySpectrum::calculate(const Forest &forest) {
  // Calculate seg_sites even if it is not a summary statistic of its own
  if (seg_sites_->position() != forest.next_base()) seg_sites_->calculate(forest);
  assert(seg_sites_->position() == forest.next_base()); 

  size_t haplotype;

  for (size_t i = at_mutation_; i < seg_sites_->countMutations(); ++i) { 
    haplotype = 0;
    for (size_t j = 0; j < seg_sites_->getHaplotype(i)->size(); ++j) 
      haplotype += (*(seg_sites_->getHaplotype(i)))[j];
    sfs_.at(haplotype - 1) += 1; 
  }

  at_mutation_ = seg_sites_->countMutations();
}

void FrequencySpectrum::printLocusOutput(std::ostream &output) const {
  output << "SFS: " << sfs_ << std::endl;
}

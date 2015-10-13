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

#include "tmrca.h"

void TMRCA::calculate(const Forest &forest) {
  if (forest.calcSegmentLength() == 0) return;
  tmrca_.push_back(forest.getTMRCA(true));
  tree_length_.push_back(forest.getLocalTreeLength(true));
}


void TMRCA::printLocusOutput(std::ostream &output) const {
  for (size_t i = 0; i < tmrca_.size(); ++i) {
    output << "time:\t" << tmrca_.at(i) << " \t" << tree_length_.at(i) << "\n";
  }
}

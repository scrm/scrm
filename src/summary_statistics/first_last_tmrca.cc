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

#include "first_last_tmrca.h"

void FirstLastTMRCA::calculate(const Forest &forest) {
  if (forest.current_base() == 0.0) first_tmrca_ = forest.getTMRCA(true);
  if (forest.next_base() == forest.model().loci_length()) last_tmrca_ = forest.getTMRCA(true); 
}

void FirstLastTMRCA::printLocusOutput(std::ostream &output) {
  output << "FLTMRCA: " << first_tmrca_ << " " << last_tmrca_ << std::endl;
}

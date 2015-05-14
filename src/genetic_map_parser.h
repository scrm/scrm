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

#ifndef scrm_src_genetic_map_parser
#define scrm_src_genetic_map_parser

#include <iostream>
#include <fstream>
#include <model.h>

class GeneticMapParser {
 public:
  GeneticMapParser() {

  }

  void add_to_model(Model &model); 

 private:
  std::ifstream map_file_;
  size_t line_offset_;
  size_t pos_col_;
  size_t ratio_col_;
}

#endif

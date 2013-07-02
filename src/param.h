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

#ifndef scrm_param
#define scrm_param

#include <iostream>
#include <iomanip>      
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>  
#include <stdio.h>
#include <stdexcept>
#include <boost/lexical_cast.hpp> 
#include "model.h"

class Param {
 public:
  friend std::ostream& operator<< (std::ostream& stream, const Param& param) {
    stream << "scrm ";
    for (size_t i = 1; i < param.argc_; ++i) {
      stream << param.argv_[i] << " ";
    }
    return stream;
  };

  int random_seed;
  std::string treefile;
  
  std::string tmrca_NAME;
  
  void init();

  bool seg_bool() const { return seg_bool_; }
  bool tmrca_bool() const { return tmrca_bool_; }

  void set_seg_bool(const bool &seg_bool) { seg_bool_ = seg_bool; } 
  void set_tmrca_bool(const bool &tmrca_bool) { tmrca_bool_ = tmrca_bool; } 

  Param() : argc_(0), argv_(NULL) { };
  Param(int argc, char *argv[]) : argc_(argc), argv_(argv) { }

  Model* parse();
  
  void print_param();
  
  bool tree_bool;
  bool log_bool;
  std::string log_NAME;
  
  void log_param();
   
 private:
  const int argc_;
  char * const* argv_;

  bool tmrca_bool_;
  bool seg_bool_; 
};

void print_help();
void print_example();
void print_option();	

template<class T>
T readInput(char input[])
{
  return boost::lexical_cast<T>(input);
}

#endif

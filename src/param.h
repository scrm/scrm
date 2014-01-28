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
 #ifdef UNITTEST
  friend class TestParam;
 #endif
  // Constructors
  Param() : argc_(0), argv_(NULL) { };
  Param(int argc, char *argv[], bool directly_called=true) : argc_(argc), argv_(argv) , directly_called_(directly_called) { }

  // Getters and setters
  bool seg_bool() const { return seg_bool_; }
  void set_seg_bool(const bool &seg_bool) { seg_bool_ = seg_bool; } 

  size_t random_seed() const { return random_seed_; }
  void set_random_seed(const int seed ){ this->random_seed_ = seed; }

  // Other methods
  void init();
  void log_param();
  friend std::ostream& operator<< (std::ostream& stream, const Param& param);

  void parse(Model &model);
  void nextArg(std::string option);
  void print_param();

  // Member variables
  
  bool tree_bool;
  bool tmrca_bool;
  bool log_bool;

  std::string tree_NAME;
  std::string tmrca_NAME;
  std::string log_NAME;

  template<class T>
  T readNextInput() {
    ++argc_i;

    if (argc_i >= argc_) {
      throw std::invalid_argument(std::string("Not enough parameters when parsing options"));
    }

    return boost::lexical_cast<T>(argv_[argc_i]);
  }

 private:
  const int argc_;
  int argc_i;
  char * const* argv_;
  bool seg_bool_;
  size_t random_seed_;
  bool directly_called_;
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

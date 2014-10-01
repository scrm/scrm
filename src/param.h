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

#ifndef scrm_param
#define scrm_param

#include <vector>
#include <iostream>
#include <iomanip>      
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>  
#include <stdio.h>
#include <stdexcept>
#include <memory>

#include "model.h"
#include "summary_statistics/summary_statistic.h"
#include "summary_statistics/tmrca.h"
#include "summary_statistics/first_last_tmrca.h"
#include "summary_statistics/seg_sites.h"
#include "summary_statistics/frequency_spectrum.h"
#include "summary_statistics/newick_tree.h"
#include "summary_statistics/oriented_forest.h"

class Param {
 public:
 #ifdef UNITTEST
  friend class TestParam;
 #endif

  // Constructors
  Param() : argc_(0), argv_(NULL) { init(); }
  Param(const std::string &arg);
  Param(int argc, char *argv[], bool directly_called=true) : 
      argc_(argc), argv_(argv), directly_called_(directly_called) {
        init();
  }

  ~Param() {
    for (char* arg : argv_vec_) delete arg;
  }
 
  /** Move Operator */
  Param(Param&& other) {
    argc_ = other.argc_;
    argc_i = other.argc_i;
    argv_ = other.argv_;
    random_seed_ = other.random_seed_;  
    directly_called_ = other.directly_called_;
    help_ = other.help_;
    version_ = other.version_;
    std::swap(argv_vec_, other.argv_vec_);
  }

  /** Copy Assignment Operator */
  Param& operator=(Param other) {
    argc_ = other.argc_;
    argc_i = other.argc_i;
    argv_ = other.argv_;
    random_seed_ = other.random_seed_;  
    directly_called_ = other.directly_called_;
    help_ = other.help_;
    version_ = other.version_;
    std::swap(argv_vec_, other.argv_vec_);
    return *this;
  }

  void init() {
    this->set_random_seed(-1);
    this->set_help(false);
    this->set_version(false);
    argc_i = 0;
  }

  // Getters and setters
  bool help() const { return help_; }
  bool version() const { return version_; }
  size_t random_seed() const { return random_seed_; }
  void set_random_seed(const size_t seed) { this->random_seed_ = seed; }

  // Other methods
  void printHelp(std::ostream& stream);

  friend std::ostream& operator<< (std::ostream& stream, const Param& param);

  void parse(Model &model);

  template<class T>
  T readNextInput() {
    ++argc_i;

    if (argc_i >= argc_) {
      throw std::invalid_argument(std::string("Not enough parameters when parsing options: ") + argv_[argc_i-1]);
    }

    T input;
    std::stringstream ss(argv_[argc_i]);
    ss >> input;
    if (ss.fail()) {
      throw std::invalid_argument(std::string("Failed to parse option: ") + argv_[argc_i]);
    }
    return input;
  }

 private:
  Param(const Param &other);
  
  void set_help(const bool help) { this->help_ = help; } 
  void set_version(const bool version) { this->version_ = version; } 

  int argc_;
  int argc_i;
  char * const* argv_;
  size_t random_seed_;  
  bool directly_called_;
  bool help_;
  bool version_;
  std::vector<char*> argv_vec_;
};
#endif

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

#ifndef scrm_param
#define scrm_param

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <random>
#include <iterator>
#include <sstream>


#include "model.h"
#include "summary_statistics/summary_statistic.h"
#include "summary_statistics/tmrca.h"
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
  Param() { init(); }

  Param(const std::string &arg) { 
    std::istringstream iss(arg);
    copy(std::istream_iterator<std::string>(iss),
         std::istream_iterator<std::string>(),
         std::back_inserter(argv_));
    directly_called_ = true;
    init();
  }

  Param(int argc, char *argv[], bool directly_called=true) : 
    directly_called_(directly_called) {
    argv_ = std::vector<std::string>(argv + 1, argv + argc);
    init();
  }

  void init() {
    this->seed_set_ = false;
    this->random_seed_ = 0;
    this->set_help(false);
    this->set_version(false);
    this->set_precision(6);
    this->set_print_model(false);
    this->argv_i = argv_.begin();
  }

  // Getters and setters
  bool help() const { return help_; }
  bool version() const { return version_; }
  bool read_init_genealogy() const { return this->read_init_genealogy_; }
  size_t random_seed() const { return random_seed_; }
  size_t precision() const { return precision_; }
  bool seed_is_set() const { return this->seed_set_; }
  bool print_model() const { return this->print_model_; }

  void set_precision ( const size_t p ) { this->precision_ = p; }
  void set_random_seed(const size_t seed) { 
    this->random_seed_ = seed;
    this->seed_set_ = true; 
  }
  void set_print_model(const bool print_model) { print_model_ = print_model; }

  // Other methods
  void printHelp(std::ostream& stream);

  friend std::ostream& operator<< (std::ostream& stream, const Param& param);

  Model parse();

  template<class T>
  T readNextInput() {
    ++argv_i;

    if (argv_i == argv_.end()) {
      throw std::invalid_argument(std::string("Unexpected end of arguments"));
    }

    return convert<T>(*argv_i);
  }

  template<class T>
  T convert(const std::string &arg) {
    T value;
    std::stringstream ss(arg);
    ss >> value;
    if (ss.fail()) {
      throw std::invalid_argument(std::string("Failed to parse option: ") + arg);
    }
    return value;
  }

  // Read to double first and then cast to int to support scientific notation
  size_t readNextInt() {
    return readNextInput<double>();
  }


  std::vector < std::string > init_genealogy;

 private:
  void set_help(const bool help) { this->help_ = help; } 
  void set_version(const bool version) { this->version_ = version; } 

  std::vector<std::string> argv_;
  std::vector<std::string>::iterator argv_i;
  size_t seed_set_;
  size_t random_seed_;
  size_t precision_;
  bool directly_called_;
  bool help_;
  bool version_;
  bool read_init_genealogy_;
  bool print_model_;
};
#endif

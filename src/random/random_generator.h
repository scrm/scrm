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

#ifndef scrm_src_random_random_generator
#define scrm_src_random_random_generator

#include "../macros.h" // Needs to be before cassert

#include <cassert>
#include <cmath>
#include <memory>

#include "fastfunc.h"


class RandomGenerator
{
 friend class MersenneTwister;
 public:
  RandomGenerator() : ff_(std::make_shared<FastFunc>()) { };
  RandomGenerator(std::shared_ptr<FastFunc> ff) : ff_(ff) { };

  virtual ~RandomGenerator() {}

  //Getters & Setters
  size_t seed() const { return seed_; }

  virtual void set_seed(const size_t seed) {
    this->seed_ = seed;
  }

  virtual double sample() =0;

  // Base class methods
  // Initialize unit_exponential; must be called when the random generator is up and running
  void initializeUnitExponential() {
    this->unit_exponential_ = sampleUnitExponential();
  }

  // Uniformly samples a number out of 0, ..., range-1
  // Unit tested
  int sampleInt(const int max_value) {
    return(static_cast<int>(this->sample()*max_value));
  }

  // Samples from an exponential distribution
  // Unit tested
  double sampleExpo(const double lambda) {
    return sampleUnitExponential() / lambda;
  }

  // Samples from an exponential distribution; return -1 if beyond limit
  // If a limit is known, this version is faster than the standard one
  // Unit tested
  double sampleExpoLimit(const double lambda, const double limit) {
    return sampleExpoExpoLimit(lambda, 0, limit);
  }

  double sampleExpoExpoLimit(const double b, const double c, const double limit);

#ifdef UNITTEST
  friend class TestRandomGenerator;
#endif

  // fast functions
  std::shared_ptr<FastFunc> ff() { return this->ff_; }


 protected:
  // Sample from a unit exponential distribution
  // Unit tested
  // fastlog can return 0, which causes a bug in scrm.
  // log or fastlog does seems to have an influence on the runtime.
  virtual double sampleUnitExponential() {
    return -(this->ff()->fastlog(sample()));
    //return -std::log(sample());
  }

 protected:
  // seed
  size_t seed_;
  // cache for a unit-exponentially distributed variable
  double unit_exponential_;
  std::shared_ptr<FastFunc> ff_;
};

#endif

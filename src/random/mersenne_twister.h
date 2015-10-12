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

#ifndef scrm_src_random_mersenne_twister
#define scrm_src_random_mersenne_twister

#include <random>
#include <memory>
#include "random_generator.h"

class MersenneTwister : public RandomGenerator
{
 public:
  MersenneTwister();
  MersenneTwister(const size_t seed);
  MersenneTwister(const bool use_seed, size_t seed);

  MersenneTwister(std::shared_ptr<FastFunc> ff):RandomGenerator(ff) {
      this->construct_common(generateRandomSeed());
  }

  MersenneTwister(const size_t seed, std::shared_ptr<FastFunc> ff):RandomGenerator(ff) {
      this->construct_common(seed);
  }

  ~MersenneTwister() {};

  void set_seed(const size_t seed);
  void construct_common(const size_t seed);

  double sample() { return unif_(mt_); }

 protected:
  std::mt19937_64 mt_; 
  std::uniform_real_distribution<> unif_;

 private:
  size_t generateRandomSeed() const;
};

#endif

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

#include "mersenne_twister.h"

void MersenneTwister::construct_common(const size_t seed){
  unif_ = std::uniform_real_distribution<>(0, 1);
  this->set_seed(seed);
}

MersenneTwister::MersenneTwister() {
  this->construct_common(generateRandomSeed());
}

MersenneTwister::MersenneTwister(const size_t seed){
  this->construct_common(seed);
}

MersenneTwister::MersenneTwister(const bool use_seed, size_t seed){
  if (!use_seed) seed = generateRandomSeed();
  this->construct_common(seed);
}


/**
 * @brief Generates a random seed using entropy provided by the operating
 * system. 
 *
 * @return A random int between 0 and 2^32
 */
size_t MersenneTwister::generateRandomSeed() const {
  std::random_device rd;
  std::uniform_int_distribution<size_t> dist(0, 4294967295); // 0 - 2^32-1
  return(dist(rd));
}

void MersenneTwister::set_seed(const size_t seed) {
  RandomGenerator::set_seed(seed);
  mt_ = std::mt19937_64(seed);
  this->initializeUnitExponential();
}


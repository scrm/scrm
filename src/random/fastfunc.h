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

#ifndef scrm_src_random_fastfunc
#define scrm_src_random_fastfunc

#include <stdio.h>

#include <stdint.h>
#include <math.h>
//#include <malloc.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <iostream>

class FastFunc {
 public:
  FastFunc();
  ~FastFunc();

  // Methods
  double fastlog(double);              /* about as fast as division; about as accurate as logf */
  double fastexp_up(double y);  /* upper bound to exp; at most 6.148% too high.  10x as fast as exp */
  double fastexp_lo(double y);  /* lower bound to exp; at most 5.792% too low.  10x as fast as exp */

 protected:
  float* build_fastlog_table();
  
  static constexpr double EXP_A=1048576/M_LN2;
  static constexpr long long EXP_C_LO=90254;
  static constexpr long long EXP_C_UP=-1;

  float* fastlog_table_;
};


class _float_bits {
 public:
  union {
    float f;
    int32_t i;
  } n;
  _float_bits(float f) { n.f=f; }
  _float_bits(int32_t i) { n.i=i; }
  float f() const { return n.f; }
  int32_t i() const { return n.i; }
};


// Inline definitions
inline FastFunc::FastFunc() {
  fastlog_table_ = build_fastlog_table();
}

inline FastFunc::~FastFunc() {
    free(fastlog_table_);
}

// Fast and fairly tight upper and lower bounds for exp(x)
// A quick test showed a ~10x speed improvement
// See: Nicol N. Schraudolf, A Fast, Compact Approximation of the Exponential Function, Neural Computation 11, 853-862 (1999)
// http://nic.schraudolph.org/pubs/Schraudolph99.pdf


inline double FastFunc::fastexp_up(double y) {
  if (y<-700) return 0.0;
  if (y>700) return INFINITY;
  union {
    double d;
    int64_t i;
  } n;
  //int64_t* yint = (int64_t*)(&y);
  //*yint = (((long long)(EXP_A*y)) + (1072693248 - EXP_C_UP)) << 32;
  n.i = (((long long)(EXP_A*y)) + (1072693248 - EXP_C_UP)) << 32;
  //return *(double*)(yint);
  return n.d;
}

inline double FastFunc::fastexp_lo(double y) {
  if (y<-700) return 0.0;
  if (y>700) return INFINITY;
  union {
    double d;
    int64_t i;
  } n;
  //int64_t* yint = (int64_t*)(&y);
  //*yint = (((long long)(EXP_A*y)) + (1072693248 - EXP_C_LO)) << 32;
  //return *(double*)(yint);
  n.i = (((long long)(EXP_A*y)) + (1072693248 - EXP_C_LO)) << 32;
  return n.d;
}

inline double FastFunc::fastlog(double x) {
  const float offset = 2047; // as int32_t: 0x44ffe000;
  float y = x;
  int32_t* yint = (int32_t*)(&y);
  int expon = ((*yint) >> 23) - 127;   // base-2 exponent of float
  int index = ((*yint) >> 13) & 1023;  // upper 10 bits of mantissa
  *yint |= 0x7fffe000;                 // convert float into remainder of mantissa; and
  *yint &= 0x44ffffff;                 // modify exponent to get into proper range
  return (expon * M_LN2 +              // contribution of base-2 log
	  fastlog_table_[index] +      // table lookup
	  (fastlog_table_[index+1] - fastlog_table_[index]) * (*(float*)(yint) - offset) );  // linear interpolation
}


#endif

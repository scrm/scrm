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

#include <iostream>
#include <cmath>

//#include <malloc.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif


// Number of interpolation points.  If this is changed, several constants in fastlog must also be changed.
#define SIZE_DOUBLE 1024

class FastFunc {
 public:
  FastFunc();
  ~FastFunc();

  // Methods
  double fastlog(double);       /* about as fast as division; about as accurate as logf */
  double fastexp_up(double y);  /* upper bound to exp; at most 6.148% too high.  10x as fast as exp */
  double fastexp_lo(double y);  /* lower bound to exp; at most 5.792% too low.  10x as fast as exp */

 protected:
  double* build_fastlog_double_table( int );
  
  static constexpr double LN2 = 0.693147180559945309417; //ln(2)
  static constexpr double EXP_A = 1048576/LN2;
  static constexpr long long EXP_C_LO = 90254;
  static constexpr long long EXP_C_UP = -1;

  double* fastlog_double_table_;
};

// Inline definitions
inline FastFunc::FastFunc() {
  fastlog_double_table_ = build_fastlog_double_table( SIZE_DOUBLE );
}

inline FastFunc::~FastFunc() {
    free(fastlog_double_table_);
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
  n.i = (((long long)(EXP_A*y)) + (1072693248 - EXP_C_UP)) << 32;
  return n.d;
}

inline double FastFunc::fastexp_lo(double y) {
  if (y<-700) return 0.0;
  if (y>700) return INFINITY;
  union {
    double d;
    int64_t i;
  } n;
  n.i = (((long long)(EXP_A*y)) + (1072693248 - EXP_C_LO)) << 32;
  return n.d;
}

inline double FastFunc::fastlog(double x) {
  const float offset = 2047;                // as int64_t: 0x409ffc00000....
  double y = x;
  int64_t* yint = (int64_t*)(&y);
  int expon = ((*yint) >> 52) - 1023;       // base-2 exponent of float
  int index = ((*yint) >> (52-10)) & 1023;  // upper 10 bits of mantissa
  *yint |= 0x7ffffc0000000000;              // convert float into remainder of mantissa; and
  *yint &= 0x409fffffffffffff;              // modify exponent to get into proper range
  return (expon * LN2 +                   // contribution of base-2 log
	  fastlog_double_table_[index] +    // table lookup, and linear interpolation
	  (fastlog_double_table_[index+1] - fastlog_double_table_[index]) * (*(double*)(yint) - offset) );
}

#endif

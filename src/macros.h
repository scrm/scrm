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

#ifndef scrm_src_macros
#define scrm_src_macros

/* Used when building an R package (RBUILD) */
#ifdef RBUILD

// Include Rcpp Headers for Rcout.
#include "Rcpp.h"

// Suppress debug output.
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && Rcpp::Rcout

// Assure that assertions are deactivated.
#ifndef NDEBUG
#define NDEBUG
#endif

#else
/* Used for normal compilation for scrm */

// Unless compiled with options NDEBUG, we will produce a debug output using 
// 'dout' instead of cout and execute (expensive) assert statements.
#ifndef NDEBUG

// Debug mode
#ifdef UNITTEST // No debug output in unittests
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#else           // Produce debug output
#define dout std::cout
#endif

#else
// Normal Mode
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

#endif

#include <limits>
#include <cmath>
// from Knuths "The art of computer programming"
inline bool areSame(const double a, const double b, 
                    const double epsilon = std::numeric_limits<double>::epsilon()) {
  return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

#endif

#ifndef scrm_src_macros
#define scrm_src_macros

/* Used when building an R package (RBUILD) */
#ifdef RBUILD

// Include Rcpp Headers for Rcout.
#include "Rcpp.h"

// Use Rcout for debug output.
#ifndef NDEBUG
#define dout Rcpp::Rcout
#else
#define dout 0 && Rcpp::Rcout
#endif

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
#define dout 0 && std::cout
#else           // Produce debug output
#define dout std::cout
#endif

#else
// Normal Mode
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

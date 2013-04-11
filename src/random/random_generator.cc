#include "random_generator.h"

#include <iostream>


// Sample from a unit exponential distribution
double RandomGenerator::sampleUnitExponential(void) {
  double exposample = -ff.fastlog( sample() );
  return exposample;
}

// Sets new seed, and initializes unit_exponential
void RandomGenerator::set_seed(const int &seed){
  this->seed_ = seed;
  unit_exponential_ = sampleUnitExponential();
}


// Samples from an exponential distribution 
// Distribution checked -Paul
double RandomGenerator::sampleExpo(double lambda){
  return sampleUnitExponential() / lambda;
}


// Samples from an exponential distribution; return -1 if beyond limit
// If a limit is known, this version is faster than the standard one
double RandomGenerator::sampleExpoLimit(double lambda, double limit){
  if (unit_exponential_ >= limit * lambda) {
    unit_exponential_ -= limit * lambda;
    return -1;
  } else {
    double result = unit_exponential_ / lambda;
    unit_exponential_ = sampleUnitExponential();
    return result;
  }
}


// Samples waiting time, with limit, for a process with an exponentially changing rate:
//  rate(t) = b exp( c t )
// This code allows c=0, and falls back to a standard exponential if so
// For a p ~ unif(0,1), a waiting time sample is (1/c) log[ 1 - (c/b) log p ]
// It returns -1 if no event occurred; this can happen even if limit == +infinity (if c<0)
double RandomGenerator::sampleExpoExpoLimit(double b, double c, double limit){
  assert (b>0);
  assert (limit>=0);
  // for any c, the no-event condition (t=maximum waiting time) is
  //  (b/c) (exp(c t)-1) < -log p
  // For c<0 and c>0 respectively this becomes
  //  -c log p < b (exp(c t)-1)
  //  -c log p > b (exp(c t)-1)  resp.
  // These are implied by the conditions
  //  -c log p < b (exp_up(c t)-1)    [c<0]
  //  -c log p > b (exp_lo(c t)-1)    [c>0]
  // where exp_up and exp_lo are upper and lower bounds for the exp(x) function resp.
  if (c < 0) {
    double c_logp_limit = b*(ff.fastexp_lo(c*limit)-1);  // negative
    if (c*unit_exponential_ < c_logp_limit) {
      unit_exponential_ -= c_logp_limit / c;
      return -1;
    }
    double y = 1.0 + c*unit_exponential_ / b;
    unit_exponential_ = sampleUnitExponential();
    if (y <= 0.0) return -1; // no event at all
    y = ff.fastlog( y )/c;       
    if (y > limit) return -1;  // the event time; can still be beyond limit
    return y;
  } else if (c > 0) {
    double c_logp_limit = b*(ff.fastexp_up( c*limit )-1);  // positive
    if (c*unit_exponential_ > c_logp_limit) {
      unit_exponential_ -= c_logp_limit / c;
      return -1;
    }
    double y = ff.fastlog( 1.0 + c*unit_exponential_ / b ) / c;
    unit_exponential_ = sampleUnitExponential();
    if (y > limit) return -1;
    return y;
  } else {
    if (unit_exponential_ >= limit * b) {
      unit_exponential_ -= limit * b;
      return -1;
    } else {
      double result = unit_exponential_ / b;
      unit_exponential_ = sampleUnitExponential();
      return result;
    }
  }
}

    
// Uniformly samples a number out of 0, ..., range-1
// Distribution checked -Paul
int RandomGenerator::sampleInt(int range) {
  return(static_cast<int>(this->sample()*range));
}


// Uniformly samples a pair of different elements of 0, ..., range-1
// Distribution checked -Paul
void RandomGenerator::sampleTwoElements(int range, int *sample1, int *sample2) {
  *sample1 = this->sampleInt(range);
  *sample2 = this->sampleInt(range-1);
  if (*sample1 == *sample2) *sample2 = range - 1;
  assert( 0 <= *sample1 && *sample1 < range );
  assert( 0 <= *sample2 && *sample2 < range );
}




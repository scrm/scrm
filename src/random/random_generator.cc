#include "random_generator.h"


void RandomGenerator::set_seed(const int &seed){
  this->seed_ = seed;
  this->unit_exponential_ = -std::log(this->sample());
}


// Samples from an exponential distribution 
// Distribution checked -Paul
double RandomGenerator::sampleExpo(double lambda){
  return -std::log(this->sample()) / lambda;
}


// Samples from an exponential distribution; return -1 if beyond limit
// If a limit is known, this version is faster than the standard one
double RandomGenerator::sampleExpoLimit(double lambda, double limit){
  if (unit_exponential >= limit * lambda) {
    unit_exponential -= limit * lambda;
    return -1;
  } else {
    double result = unit_exponential / lambda;
    unit_exponential = -std::log(this->sample());
    return result;
  }
}


// Samples waiting time, with limit, for a process with an exponentially changing rate:
//  rate(t) = b exp( c t )
// For a p ~ unif(0,1), a waiting time sample is (1/c) log[ 1 - (c/b) log p ]
// It returns -1 if no event occurred; this can happen even if limit == +infinity (if c<0)
double RandomGenerator::sampleExpoExpoLimit(double b, double c, double limit){
  assert (b>0);
  assert (c!=0);
  // for any c, the no-event condition (t=maximum waiting time) is
  //  (b/c) (exp(c t)-1) < -log p
  // For c<0 and c>0 respectively this becomes
  //  -c log p < b (exp(c t)-1)
  //  -c log p > b (exp(c t)-1)  resp.
  // These are implied by the conditions
  //  -c log p < b (exp_up(c t)-1)    [c<0]
  //  -c log p > b (exp_lo(c t)-1)    [c>0]
  // where exp_up and exp_lo are upper and lower bounds for the exp(x) function resp.
  if (c<0) {
    double c_logp_limit = b*(EXP_UP(c*limit)-1);  // negative
    if (c*unit_exponential < c_logp_limit) {
      unit_exponential -= c_logp_limit / c;
      return -1;
    }
    double y = 1.0 + c*unit_exponential / b;
    unit_exponential = -std::log(this->sample());
    if (y <= 0.0) return -1; // no event at all
    y = std::log(y)/c;       
    if (y>limit) return -1;  // the event time; can still be beyond limit
    return y;
  } else {
    double c_logp_limit = b*(EXP_UP(c*limit)-1);  // positive
    if (c*unit_exponential > c_logp_limit) {
      unit_exponential -= c_logp_limit / c;
      return -1;
    }
    double y = std::log( 1.0 + c*unit_exponential / b ) / c;
    unit_exponential = -std::log(this->sample());
    if (y>limit) return -1;
    return y;
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



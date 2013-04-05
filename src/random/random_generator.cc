#include "random_generator.h"

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


// Samples from an exponential distribution 
// Distribution checked -Paul
double RandomGenerator::sampleExpo(double lambda){
  return -std::log(this->sample()) / lambda;
}

void RandomGenerator::set_seed(const int &seed){
  this->seed_ = seed;
}

#include "mersenne_twister.h"

MersenneTwister::MersenneTwister() {
  this->set_seed(std::time(0));
};

MersenneTwister::MersenneTwister(int seed){
  this->set_seed(seed);
}

MersenneTwister::~MersenneTwister() { } ;
                         
double MersenneTwister::sample() {
  return(unif(rng));
}

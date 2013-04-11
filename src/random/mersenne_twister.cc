#include "mersenne_twister.h"

MersenneTwister::MersenneTwister() {
  this->set_seed(std::time(0));
  //std::cout << "Seed: " << this->seed() << std::endl;
};

MersenneTwister::MersenneTwister(int seed){
  this->set_seed(seed);
  //std::cout << "Seed: " << this->seed() << std::endl;
}

MersenneTwister::~MersenneTwister() { } ;

double MersenneTwister::sample() {
  return(unif(rng));
}

void MersenneTwister::set_seed(const int &seed) {
  this->seed_ = seed;
  this->rng.seed(static_cast<unsigned int>(seed));
}


//////////////////////////////////////////////////////////////////////////////////////////////////////

//#define UNIT_TEST
#ifdef UNIT_TEST

#include <stdio.h>

int main(int argc, char** argv) {
  
  // Sample from exponential and other distributions, and use the Kolmogorov-Smirnov test to
  // see if the correct distribution is obtained.

  const int CASES=10;
  const int SAMPLES=10000000;
  const double rate[CASES] =   {1.0,      1.0, 1.0, 1.0,      1.0, 1.0, 10.0, 1.0,      1.0,  1.0};
  const double growth[CASES] = {0.0,      0.0, 0.0, 1.0,      1.0, 1.0, 10.0, -1.0,     -1.0, -1.0};
  const double limit[CASES] =  {INFINITY, 1.0, 0.1, INFINITY, 1.0, 0.1, 0.1,  INFINITY, 1.0,  0.1};

  class MersenneTwister rg;
  double* samples = (double*)malloc( SAMPLES * sizeof(double) );
  int i,j,expire_count;
  double s;

  assert(samples);

  printf("Rate\t\tGrowth\t\tLimit\t\tExpire_z\tZtest\tKSstat\t\tKStest\n");

  for (i=0; i<CASES; i++) {
    expire_count = 0;
    for (j=0; j<SAMPLES; j++) {
      s = rg.sampleExpoExpoLimit( rate[i], growth[i], limit[i] );
      if (s==-1) {
	expire_count += 1;
	j--;
      } else {
	samples[j] = s;
      }
    }
    double expire_probability;
    if (growth[i] != 0.0) {
      expire_probability = exp( -(rate[i]/growth[i])*(exp(growth[i]*limit[i]) - 1) );
    } else {
      expire_probability = exp( -rate[i]*limit[i] );
    }
    double expire_zscore = (expire_count - (SAMPLES+expire_count)*expire_probability) / sqrt( (SAMPLES+expire_count)*expire_probability*(1.0-expire_probability) );
    std::sort( samples, samples+SAMPLES );
    double kolmogorov_smirnov_d = 0.0;
    for (j=0; j<SAMPLES; j++) {
      double x = samples[j];
      double empirical_lower = (double)j/SAMPLES;
      double empirical_upper = (double)(j+1)/SAMPLES;
      double cdf;
      if (growth[i] != 0.0) {
	cdf = 1.0 - (exp( -(rate[i]/growth[i])*(exp(growth[i]*x) - 1.0) ) - expire_probability) / (1.0 - expire_probability);
      } else {
	cdf = 1.0 - (exp( -rate[i]*x ) - expire_probability) / (1.0 - expire_probability);
      }
      double deviation = std::max( std::abs( cdf - empirical_lower ), std::abs( cdf - empirical_upper ) );
      kolmogorov_smirnov_d = std::max( kolmogorov_smirnov_d, deviation );
    }
    kolmogorov_smirnov_d *= sqrt(SAMPLES);
    double ks_critical_value = 1.63;  // 1% level
    double z_critical_value = 2.58;   // 1% level, two-sided
    printf("%9.7e\t%9.7e\t%9.7e\t%9.7e\t%s\t%9.7e\t%s\n", 
	   rate[i], growth[i], limit[i], 
	   expire_zscore, (abs(expire_zscore)<z_critical_value)?"ok":"FAIL",
	   kolmogorov_smirnov_d, (kolmogorov_smirnov_d<ks_critical_value)?"ok":"FAIL");
  }
  return 0;
}

#endif

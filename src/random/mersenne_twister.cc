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


// If unit testing, compile as follows:
//
// g++ -c -O3 random_generator.cc mersenne_twister.cc
// g++ -O3 mersenne_twister.o random_generator.o ../icsilog/icsilog.o -lm
//
//#define UNIT_TEST
#ifdef UNIT_TEST

#include <stdio.h>

double speedtest(int oper, double rate, double growth, double limit) {

  int i;
  double x = 1.0;
  double y = 0.999;
  double z = 0.0;
  class MersenneTwister rg;
  for (i=0; i<100000000; i++) {
    x += 1e-8;
    y = x;
    switch(oper) {
    case 0: break;
    case 1: y = x+y; break;
    case 2: y = x/x; break;
    case 3: y = x*x; break;
    case 4: y = log(x); break;
    case 5: y = rg.mylog(x); break;
    case 6: y = exp(x); break;
    case 7: y = EXP_LO(x); break;
    case 8: y = rg.sample(); break;
    case 9: y = rg.sampleExpo( 1.0 ); break;
    case 10: y = rg.sampleExpoLimit( 1.0, 0.1 ); break;
    default: y = rg.sampleExpoExpoLimit( rate, growth, limit ); break;
    }
    z += y;
  }
  return z;  // to avoid optimizing everything away
}

void exptest() {

  class MersenneTwister rg;
  int i;
  for (i=0; i<10; i++) {
    double x = (rg.sample()-0.5) * 1400.0;
    double true_exp = exp(x);
    double lower_bound = EXP_LO(x);
    double upper_bound = EXP_UP(x);
    assert (lower_bound > true_exp * (1.0 - 0.05792));
    assert (upper_bound >= true_exp);
    assert (upper_bound < true_exp * (1.0 + 0.06148));
  }
}
    
int main(int argc, char** argv) {
  
  const int CASES=10;
  const int SAMPLES=5000000;
  const double rate[CASES] =   {1.0,      1.0, 1.0, 1.0,      1.0, 1.0, 10.0, 1.0,      1.0,  1.0};
  const double growth[CASES] = {0.0,      0.0, 0.0, 1.0,      1.0, 1.0, 10.0, -1.0,     -1.0, -1.0};
  const double limit[CASES] =  {INFINITY, 1.0, 0.1, INFINITY, 1.0, 0.1, 0.1,  INFINITY, 1.0,  0.1};

  class MersenneTwister rg;
  double* samples = (double*)malloc( SAMPLES * sizeof(double) );
  int i,j,expire_count;
  double s;

  assert(samples);

  if (argc == 2) {
    // with any argument, run the timing tests
    const int LLTEST=11;
    const char* testnames[LLTEST+1] = {"nop\t","+\t","/\t","*\t","log()\t","mylog()\t","exp()\t","fast_exp()","raw sample","sample(1.0)","sample(1.0,0.1)","sample\t"};
    printf("Test\t\tRate\tGrowth\tLimit\tTime\n");
    for (int i=0; i<LLTEST + CASES; i++) {
      clock_t start = clock(), diff;
      double r = i>=LLTEST ? rate[i-LLTEST] : 0, g = i>=LLTEST ? growth[i-LLTEST] : 0, l = i>=LLTEST ? limit[i-LLTEST] : 0;
      speedtest(i, r, g, l);
      diff = clock() - start;
      printf("%s\t%1.4f\t%1.4f\t%1.4f\t%ld\n",testnames[i>=LLTEST ? LLTEST : i],r,g,l,diff * 1000 / CLOCKS_PER_SEC);
    }
    return 0;
  }
      
  exptest();

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

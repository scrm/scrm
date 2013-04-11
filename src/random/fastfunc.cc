#include "fastfunc.h"
#include "mersenne_twister.h"

float* FastFunc::build_fastlog_table() {
  float * table = (float*) malloc( 1025*sizeof(float) );
  for (int index=0; index<1025; index++) {
    double dlog1 = log( 1.0 + (index) / 1024. );
    double dlog2 = log( 1.0 + (index+1.0) / 1024. );
    double middiff1 = log( 1.0 + (index+0.5)/1024. ) - (dlog1 + (dlog2-dlog1)*0.5);
    table[index] = (float)(dlog1 + middiff1/2);
  }
  return table;
}



// If unit testing, compile as follows:
//
// g++ -c -O3 random_generator.cc mersenne_twister.cc
// g++ -O3 mersenne_twister.o random_generator.o ../icsilog/icsilog.o -lm
//
//#define UNIT_TEST
#ifdef UNIT_TEST

#include <stdio.h>

void logtest() {

  class MersenneTwister rg;
  class FastFunc ff;
  double d = 1.0;
  double maxdiff = 0.0;

  for (int i=0; i<10000000; i++) {
    double logd = ff.fastlog( d );
    double truelogd = log( d );
    maxdiff = std::max( maxdiff, std::abs( logd-truelogd ) );
    d += 1e-7;
  }
  assert( maxdiff < 1.2e-7 );

  printf("fastlog OK\n");

}


void exptest() {

  class MersenneTwister rg;
  class FastFunc ff;
  int i;

  for (i=0; i<10000000; i++) {
    double x = (rg.sample()-0.5) * 1400.0;
    double true_exp = exp(x);
    double lower_bound = ff.fastexp_lo(x);
    double upper_bound = ff.fastexp_up(x);
    //std::cout << x << " " << true_exp << " " << lower_bound << " " << upper_bound << std::endl;
    assert (lower_bound <= true_exp);
    assert (lower_bound > true_exp * (1.0 - 0.05792));
    assert (upper_bound >= true_exp);
    assert (upper_bound < true_exp * (1.0 + 0.06148));
  }
  printf("fastexp_lo/hi OK\n");
}

int main(int argc, char **argv) {
  logtest();
  exptest();
}
    
#endif


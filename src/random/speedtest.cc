#include "fastfunc.h"
#include "mersenne_twister.h"


double speedtest(int oper, double rate, double growth, double limit) {

  int i;
  double x = 1.0;
  float xf = 1.0;
  double y = 0.999;
  float yf = 0.0;
  double z = 0.0;
  float zf = 0.0;

  class MersenneTwister rg;
  class FastFunc ff;

  for (i=0; i<100000000; i++) {
    x += 1e-8;
    xf = float(x);
    y = x;
    switch(oper) {
    case 0: break;
    case 1: y = x+y; break;
    case 2: y = x/x; break;
    case 3: y = x*x; break;
    case 4: y = log(x); break;
    case 5: yf = logf(xf); break;
    case 6: y = ff.fastlog(x); break;
    case 7: y = exp(x); break;
    case 8: y = ff.fastexp_lo(x); break;
    case 9: y = rg.sample(); break;
    case 10: y = rg.sampleExpo( 1.0 ); break;
    case 11: y = rg.sampleExpoLimit( 1.0, 0.1 ); break;
    default: y = rg.sampleExpoExpoLimit( rate, growth, limit ); break;
    }
    z += y;
    zf += yf;
  }
  return z + zf;  // to avoid optimizing everything away
}



int main(int argc, char** argv) {
  
  const int CASES=10;
  const int LLTEST=12;
  const char* testnames[LLTEST+1] = {"nop\t","+\t","/\t","*\t","log()\t","logf()\t","fastlog()","exp()\t","fastexp()","raw sample","sample(1.0)","sample(1.0,0.1)","sample\t"};
  const double rate[CASES] =   {1.0,      1.0, 1.0, 1.0,      1.0, 1.0, 10.0, 1.0,      1.0,  1.0};
  const double growth[CASES] = {0.0,      0.0, 0.0, 1.0,      1.0, 1.0, 10.0, -1.0,     -1.0, -1.0};
  const double limit[CASES] =  {INFINITY, 1.0, 0.1, INFINITY, 1.0, 0.1, 0.1,  INFINITY, 1.0,  0.1};

  class MersenneTwister rg;

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


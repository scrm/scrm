#ifndef scrm_src_random_random_generator
#define scrm_src_random_random_generator

#include <cassert>
#include <cmath>
#include <stdint.h>

// remove if normal logarithm is to be used
#define ICSILOG_PRECISION 12

#ifdef ICSILOG_PRECISION
#include "../icsilog/icsilog.h"
#include <malloc.h>
#endif


class RandomGenerator
{
  public:
   RandomGenerator() { 
     unit_exponential_ = 0.5; 
#ifdef ICSILOG_PRECISION
     icsi_table_ = (float*) malloc(((int) pow(2,ICSILOG_PRECISION))*sizeof(float));
     fill_icsi_log_table2(ICSILOG_PRECISION,icsi_table_);     
#endif
   }
   virtual ~RandomGenerator() {}
 
   //Getters & Setters
   int seed() { return this->seed_; }

   //Virtual methods
   virtual void initialize() =0;
   virtual double sample() =0;
   virtual void set_seed(const int &seed);

   //Base class methods
   int sampleInt(int max_value);
   void sampleTwoElements(int size, int *sample1, int *sample2);
   double mylog(double);
   double sampleExpo(double lambda);
   double sampleExpoLimit(double lambda, double limit);
   double sampleExpoExpoLimit(double b, double c, double limit);

#ifdef UNITTEST
   friend class TestRandomGenerator;
#endif

  protected:
   // seed
   int seed_;
   // cache for a unit-exponentially distributed variable
   double unit_exponential_;

#ifdef ICSILOG_PRECISION
   float *icsi_table_;
#endif

};


// Fast and pretty tight upper and lower bounds for exp(x)
// A quick test showed a ~10x speed improvement
// See: Nicol N. Schraudolf, A Fast, Compact Approximation of the Exponential Function, Neural Computation 11, 853-862 (1999)
// http://nic.schraudolph.org/pubs/Schraudolph99.pdf

class _DoubleBits {
  union {
    double d;
    struct {
      int32_t j;
      int32_t i;
    } n;
  } _double_bits;
public:
  _DoubleBits( int i ) {
    _double_bits.n.i = i;
    _double_bits.n.j = 0;
  }
  double value() const { return _double_bits.d; }
};

const double EXP_A=1048576/M_LN2;
const int EXP_C_LO=90253;  /* for lower bound to exp */
const int EXP_C_UP=-1;     /* for upper bound to exp */

// upper bound to exp; at most 6.148% too high
inline double EXP_UP(double y) {
  return (y<-700.0 ? 0.0 : (y>700.0 ? INFINITY : _DoubleBits(EXP_A*y + (1072693248 - EXP_C_UP)).value() ) ); 
}

// lower bound to exp; at most 5.792% too low
inline double EXP_LO(double y) {
  return (y<-700.0 ? 0.0 : (y>700.0 ? INFINITY : _DoubleBits(EXP_A*y + (1072693248 - EXP_C_LO)).value() ) );
}

// fast log
inline double RandomGenerator::mylog(double x) {
#ifdef ICSILOG_PRECISION
  return icsi_log_v2(x, this->icsi_table_, ICSILOG_PRECISION);
#else
  // return std::log(x);
#endif
}





#endif

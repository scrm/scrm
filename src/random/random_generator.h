#ifndef scrm_src_random_random_generator
#define scrm_src_random_random_generator

#include <cassert>
#include <cmath>

class RandomGenerator
{
  public:
   RandomGenerator() { unit_exponential_ = 0.5; }
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
   double sampleExpo(double lambda);

#ifdef UNITTEST
   friend class TestRandomGenerator;
#endif

  protected:
   // seed
   int seed_;
   // cache for a unit-exponentially distributed variable
   double unit_exponential_;
};


// Fast and pretty tight upper and lower bounds for exp(x)
// A quick test showed a ~10x speed improvement
// See: Nicol N. Schraudolf, A Fast, Compact Approximation of the Exponential Function, Neural Computation 11, 853-862 (1999)
// http://nic.schraudolph.org/pubs/Schraudolph99.pdf

class _DoubleBits {
  union {
    double d;
    struct {
#ifdef BIG_ENDIAN
    int j, i;
#else
    int i, j;
#endif
    } n;
  } _double_bits;
public:
  _DoubleBits( int i ) {
    _double_bits.n.i = i;
    _double_bits.n.j = 0;
  }
  double value() const { return _double_bits.d; }
};

#define EXP_A (1048576/M_LN2) 
#define EXP_C_UP 90253  /* for upper bound to exp */
#define EXP_C_LO -1     /* for lower bound to exp */
#define EXP_UP(y) ( _DoubleBits(EXP_A*(y) + (1072693248 - EXP_C_UP)).value() ) // upper bound to exp; at most 6.148% too high
#define EXP_LO(y) ( _DoubleBits(EXP_A*(y) + (1072693248 - EXP_C_LO)).value() ) // lower bound to exp; at most 5.792% too low


#endif

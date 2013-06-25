#ifndef scrm_src_random_random_generator
#define scrm_src_random_random_generator

#include <cassert>
#include <cmath>
#include <boost/math/distributions/poisson.hpp>

#include "fastfunc.h"


class RandomGenerator
{
 public:
  RandomGenerator() {};
  virtual ~RandomGenerator() {}
 
   //Getters & Setters
   int seed() { return this->seed_; }

   //Virtual methods
   virtual void initialize() =0;
   virtual double sample() =0;
   virtual void set_seed(const int&);

   //Base class methods
   int sampleInt(int max_value);
   void sampleTwoElements(int size, int *sample1, int *sample2);
   double sampleExpo(double lambda);
   double sampleExpoLimit(double lambda, double limit);
   double sampleExpoExpoLimit(double b, double c, double limit);

   int samplePoisson(double lambda) {
     boost::math::poisson poisson(lambda);
     return quantile(poisson, sample());
   }

#ifdef UNITTEST
   friend class TestRandomGenerator;
#endif

   //Protected methods
  protected:
   double sampleUnitExponential();

  protected:
   // seed
   int seed_;
   // cache for a unit-exponentially distributed variable
   double unit_exponential_;
   // fast functions
   FastFunc ff;

};

#endif

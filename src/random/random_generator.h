#ifndef scrm_src_random_random_generator
#define scrm_src_random_random_generator

#include <cmath>

class RandomGenerator
{
  public:
   RandomGenerator() {};
   ~RandomGenerator() {};
 
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
   int seed_;
};

#endif

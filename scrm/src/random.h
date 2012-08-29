#ifndef scrm_src_random
#define scrm_src_random

#include <vector>
#include <boost/random.hpp>
#include <ctime> 

class RandomGenerator
{
  public:
   RandomGenerator();
   RandomGenerator(int seed);
   ~RandomGenerator();
                         
   //Getters & Setters
   int seed() { return this->seed_; }

   void initialize();
   int sampleInt(int max_value);

   void sampleTwoElements(int size, int *sample1, int *sample2);
   double sample();
   double sampleExpo(double lambda);

#ifdef UNITTEST
    friend class TestRandomGenerator;
#endif

  private:
   void set_seed(const int &seed);
   int seed_;
     
   typedef boost::mt19937 rng_type;
   boost::uniform_01<double> unif;
   rng_type rng;

};

#endif

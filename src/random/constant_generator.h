#ifndef scrm_src_random_constant_generator
#define scrm_src_random_constant_generator

#include <iostream>
#include <fstream>

#include "random_generator.h"

using namespace std;

class ConstantGenerator : public RandomGenerator
{
  public:
   ConstantGenerator();
   ConstantGenerator(int seed);
   ~ConstantGenerator();

   virtual double sample();
   void initialize();
};

#endif

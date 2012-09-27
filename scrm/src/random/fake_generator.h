#ifndef scrm_src_random_fake_generator
#define scrm_src_random_fake_generator

#include <iostream>
#include <fstream>

#include "random_generator.h"

using namespace std;

class FakeRandomGenerator : public RandomGenerator
{
  public:
   FakeRandomGenerator();
   FakeRandomGenerator(int seed);
   ~FakeRandomGenerator();

   double sample();
   void initialize();
   
  private:
    ifstream* rnd_file_;
    void set_rnd_file(ifstream* rnd_file) { this->rnd_file_ = rnd_file; }
    ifstream* rnd_file() { return this->rnd_file_; }
};

#endif

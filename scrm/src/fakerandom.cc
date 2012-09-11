#include "fakerandom.h"

#include <sstream>

FakeRandomGenerator::FakeRandomGenerator(){
  this->set_seed(0);
  this->initialize();
}

FakeRandomGenerator::FakeRandomGenerator(int seed){
  this->set_seed(seed);
  this->initialize();
}

FakeRandomGenerator::~FakeRandomGenerator(){
  this->rnd_file()->close();
  delete this->rnd_file();
}

double FakeRandomGenerator::sample() {
  if (!rnd_file()->good()) {
    this->initialize();
  }
  
  string sample;
  getline(*this->rnd_file(), sample);

  std::istringstream iss(sample); 
  double d;
  iss >> d;

  return(d);
}

void FakeRandomGenerator::initialize() {
  ifstream* rnd_file = new ifstream("random.numbers");

  if (!rnd_file->good()) {
    cout << "Unable to open file 'random.numbers'" << endl;;
    exit(1);
  }
  
  for (int i=0; i < this->seed(); i++) {
    string sample;
    getline(*rnd_file, sample);
  }

  this->set_rnd_file(rnd_file);
}

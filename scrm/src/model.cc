#include "model.h"

Model::Model() {
  this->set_sample_size(0);
  this->set_population_size(10000);
  this->set_mutation_rate(10);
  this->set_recombination_rate(20);
}

Model::Model(int sample_size) {
  this->set_sample_size(sample_size);
}

Model::~Model() { };

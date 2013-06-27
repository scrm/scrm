#include "model.h"

Model::Model() {
  this->set_sample_size(0);
  this->set_population_size(10000);
  this->set_mutation_rate(10);
  this->set_recombination_rate(20);
  this->set_smc_model(false);
  this->set_exact_window_length(0);
  this->set_prune_interval(10);
}

Model::Model(int sample_size) {
  this->set_sample_size(sample_size);
  this->set_population_size(10000);
  this->set_mutation_rate(0.00001);
  this->set_recombination_rate(0.00002);
  this->set_smc_model(false);
  this->set_exact_window_length(0);
  this->set_prune_interval(10);
}

Model::Model(param user_input) {
  this->set_sample_size(user_input.nsam);
  this->set_population_size(user_input.npop);
  this->set_mutation_rate(user_input.mutation_rate_persite);
  this->set_recombination_rate(user_input.recomb_rate_persite);
  this->set_exact_window_length(user_input.exact_window_length);
  this->set_smc_model(false);
  this->set_prune_interval(10);
}

Model::~Model() { };

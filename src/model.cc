#include "model.h"

Model::Model() { 
  current_pars_ = NULL;
}


Model::Model(size_t sample_size) {
  TimeFramePars tfp = {sample_size, 10000, 0.00001, 0.00001};
  addTimeFrame(0, tfp);
  setTime(0);

  this->set_exact_window_length(0);
  this->set_prune_interval(10);
}

Model::Model(param user_input) {
/*  this->set_sample_size(user_input.nsam);
  this->set_population_size(user_input.npop);
  this->set_mutation_rate(user_input.theta);
  this->set_recombination_rate(user_input.rho);
  this->set_exact_window_length(user_input.exact_window_length);
  this->set_smc_model(false);
  this->set_prune_interval(10); */
}

Model::~Model() { };
  
void Model::addTimeFrame(const double &time, const TimeFramePars &tfp) {
  time_frames_.insert( std::pair<double, TimeFramePars>(time, tfp) );
}


void Model::setTime(const double &time) {
  current_pars_ = &(time_frames_[time]); 
}

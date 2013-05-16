#include "model.h"

Model::Model() { 
  this->set_exact_window_length(0);
  this->set_prune_interval(0);
}


Model::Model(size_t sample_size) {
  //TimeFramePars tfp = {sample_size, 10000, 0.00001, 0.00001};

  this->set_exact_window_length(0);
  this->set_prune_interval(10);
}


Model::Model(param user_input) {
  //TimeFramePars tfp = {user_input.nsam, user_input.npop, user_input.theta, user_input.rho};

  this->set_exact_window_length(user_input.exact_window_length);
  //this->set_smc_model(false);
  this->set_prune_interval(10);
}


Model::~Model() { };

/** 
 * Function to add a new change time to the model.
 *
 * It preserves the relation between the times and the *param*_list_ containers.
 * If the same time is added multiple times, it is just added once to the model,
 * but this should not make a difference when using this function.
 *
 * @param time The time that is added
 * @returns The position the time has now in the vector
 */
size_t Model::addChangeTime(double time) {
  size_t position = 0;
  if ( change_times_.size() == 0 ) {
    change_times_.push_back(time);
    sample_sizes_list_.push_back(NULL);
    population_sizes_list_.push_back(NULL);
    return position;
  }

  std::vector<double>::iterator ti; 
  for (ti = change_times_.begin(); ti != change_times_.end(); ++ti) {
    if ( *ti == time ) return position;
    if ( *ti > time ) break; 
    ++position;
  }

  change_times_.insert(ti, time);
  
  // Add Null at the right position in all parameter vectors 
  sample_sizes_list_.insert(sample_sizes_list_.begin() + position, NULL);
  population_sizes_list_.insert(population_sizes_list_.begin() + position, NULL);
  return position;
}

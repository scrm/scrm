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
  TimeFramePars tfp = {user_input.nsam, user_input.npop, user_input.theta, user_input.rho};
  addTimeFrame(0, tfp);
  setTime(0);

  this->set_exact_window_length(user_input.exact_window_length);
  //this->set_smc_model(false);
  this->set_prune_interval(10);
}


Model::~Model() { };
  

void Model::addTimeFrame(const double &time, const TimeFramePars &tfp) {
  time_frames_.insert( std::pair<double, TimeFramePars>(time, tfp) );
  time_frame_changes_.insert(time);
}


void Model::setTime(const double &time) {
  current_pars_ = &(time_frames_[time]); 
}

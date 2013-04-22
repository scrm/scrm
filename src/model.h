/*
 * model.h
 *
 * This file contains the class Model, which is a simple container object for
 * model parameters.
 *
 */

#ifndef scrm_src_model
#define scrm_src_model

#include <cstddef>
#include <map>
#include <vector>

#include "param.h"

struct TimeFramePars {
  size_t sample_size;
  size_t population_size;
  double mutation_rate;
  double recombination_rate;
};

class Model
{
  public:

#ifdef UNITTEST
  friend class TestModel;
#endif

   Model();
   Model(size_t sample_size);
   Model(param user_input);
   
   ~Model();
                         
   // Getters
   size_t sample_size() const { return current_pars_->sample_size; }
   size_t population_size() const { return current_pars_->population_size; }
   double mutation_rate() const { return current_pars_->mutation_rate; }
   double recombination_rate() const { return current_pars_->recombination_rate; }

   bool   is_smc_model() const { return this->smc_model_; }
   size_t exact_window_length() const { return exact_window_length_; }
   size_t prune_interval() const { return prune_interval_; }
  
   void set_exact_window_length(const size_t &ewl) { exact_window_length_ = ewl; }
   void set_prune_interval(const size_t &pi) { prune_interval_ = pi; }
   // 
   void setTime(const double &time);
   const std::vector<double> &change_times() { return time_frame_changes_; }; 

  private:
   std::map<double, TimeFramePars> time_frames_;
   std::vector<double> time_frame_changes_;
   TimeFramePars const *current_pars_;

   size_t exact_window_length_;
   size_t prune_interval_;
   bool smc_model_;

   void addTimeFrame(const double &time, const TimeFramePars &tfp); 
};

#endif

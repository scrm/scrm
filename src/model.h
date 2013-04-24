/*
 * model.h
 *
 * This file contains the class Model, which is a simple container object for
 * model parameters.
 *
 */
#include"param.h"
using namespace scrm; 

#ifndef scrm_src_model
#define scrm_src_model

#include <cstddef>

class Model
{
  public:
   Model();
   Model(int sample_size);
   Model(param user_input);
   
   ~Model();
                         
   // Getters
   int sample_size() const { return this->sample_size_; }
   int population_size() const { return this->population_size_; }
   double mutation_rate() const { return this->mutation_rate_; }
   double recombination_rate() const { return this->recombination_rate_; }
   bool is_smc_model() const { return this->smc_model_; }
   size_t exact_window_length() const { return exact_window_length_; }
   size_t prune_interval() const { return prune_interval_; }
  
   // Setters 
   void set_sample_size(const int &sample_size) { this->sample_size_ = sample_size; }

   void set_population_size(const int &population_size) { 
     this->population_size_ = population_size; }

   void set_mutation_rate(const double &mutation_rate) { 
     this->mutation_rate_ = mutation_rate; }

   void set_recombination_rate(const double &recombination_rate) {
     this->recombination_rate_ = recombination_rate; }

   void set_smc_model(const bool &smc_model) {
     this->smc_model_ = smc_model; }

   void set_exact_window_length(const size_t &exact_window_length) {
     this->exact_window_length_ = exact_window_length; }
   
   void set_prune_interval(const size_t &prune_interval) {
     this->prune_interval_ = prune_interval;
   }

  private:
   int sample_size_;
   int population_size_;
   double mutation_rate_;
   double recombination_rate_;

   size_t exact_window_length_;
   size_t prune_interval_;

   bool smc_model_;
};

#endif

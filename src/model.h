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
#include <vector>
#include <cfloat>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cassert>

class Param;

class Model
{
  public:

#ifdef UNITTEST
  friend class TestModel;
  friend class TestTimeInterval;
#endif
  
  friend class Param;

   Model();
   Model(size_t sample_size);

   void init();
   
   ~Model();
   
   // Default values;
   const static size_t default_pop_size = 10000;
   const static double default_growth_rate = 0.0;
   const static double default_mig_rate = 0.0;

   // Getters
   double mutation_rate() const { return mutation_rate_; }
   double recombination_rate() const { return recombination_rate_; }
   size_t loci_length() const { return loci_length_; }

   void set_mutation_rate(double rate) { mutation_rate_ = rate; }
   void set_recombination_rate(double rate) { recombination_rate_ = rate; }
   void set_loci_length(size_t length) { loci_length_ = length; }

   double growth_rate(size_t pop = 0) const {
     if (current_growth_rates_ == NULL) return default_growth_rate;
     return current_growth_rates_->at(pop);
   }
   
   size_t population_size(size_t pop = 0) const { 
     if (current_pop_sizes_ == NULL) return default_pop_size;
     return current_pop_sizes_->at(pop);
   }

   size_t sample_size() const { return sample_times_.size(); };
   size_t sample_population(size_t sample_id) const { return sample_populations_.at(sample_id); };
   double sample_time(size_t sample_id) const { return sample_times_.at(sample_id); };

   size_t exact_window_length() const { return exact_window_length_; }
   size_t prune_interval() const { return prune_interval_; }
   size_t population_number() const { return pop_number_; }
  
   
   double getCurrentTime() const { return change_times_.at(current_time_idx_); }
   double getNextTime() const { 
    if ( current_time_idx_ + 1 >= change_times_.size() ) return FLT_MAX;
    else return change_times_.at(current_time_idx_ + 1);
   }

   void set_exact_window_length(const size_t &ewl) { exact_window_length_ = ewl; }
   void set_prune_interval(const size_t &pi) { prune_interval_ = pi; }
   void set_population_number(const size_t &pop_number) { pop_number_ = pop_number; }

   void resetTime() { 
     current_pop_sizes_ = pop_sizes_list_.at(0);
     current_growth_rates_ = growth_rates_list_.at(0);
     current_time_idx_ = 0;
   };

   void increaseTime() { 
     if ( current_time_idx_ == change_times_.size() - 1) throw std::out_of_range("Model change times out of range");
     ++current_time_idx_;

     if ( pop_sizes_list_.at(current_time_idx_) != NULL ) 
       current_pop_sizes_ = pop_sizes_list_.at(current_time_idx_);
     if ( growth_rates_list_.at(current_time_idx_) != NULL ) 
       current_growth_rates_ = growth_rates_list_.at(current_time_idx_); 
   };
  
   void print(std::ostream &os) const;

   template <typename T>
   void printVector(const std::vector<T> vec, std::ostream &os) const {
     typename std::vector<T>::const_iterator it;
     for (it = vec.begin(); it != vec.end(); ++it) os << *it << " "; 
   }

   size_t loci_number() { return loci_number_; };
   void set_loci_number(size_t loci_number) { loci_number_ = loci_number; }; 
   
   


   //const std::vector<double> &change_times() const { return model_change_times_; }

  private:
   size_t addChangeTime(double time);
   void addSampleSizes(double time, const std::vector<size_t> &samples_sizes);
   void addPopulationSizes(double time, const std::vector<size_t> &population_sizes);
   void addGrowthRates(double time, const std::vector<double> &growth_rates);
   
   template <typename T>
   void deleteParList(std::vector<T*> &parList) {
    typename std::vector<T*>::iterator it;
    for (it = parList.begin(); it != parList.end(); ++it) {
      if (*it != NULL) delete *it;
    }
    parList.clear();
   }

   std::vector<double> change_times_;

   std::vector<std::vector<size_t>*> pop_sizes_list_;
   std::vector<std::vector<double>*> growth_rates_list_;
   
   std::vector<size_t> sample_populations_;
   std::vector<double> sample_times_;

   std::vector<size_t>* current_pop_sizes_;
   std::vector<double>* current_growth_rates_;
   size_t current_time_idx_;

   double mutation_rate_;
   double recombination_rate_;
   size_t pop_number_;

   size_t loci_number_;
   size_t loci_length_;

   size_t exact_window_length_;
   size_t prune_interval_;
};

#endif

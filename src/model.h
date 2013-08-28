/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013 Paul R. Staab, Sha (Joe) Zhu and Gerton Lunter
 * 
 * This file is part of scrm.
 * 
 * scrm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

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
#include <cmath>

class Param;

class Model
{
  public:

#ifdef UNITTEST
  friend class TestModel;
  friend class TestTimeInterval;
#endif
  
  friend class Param;
  friend std::ostream& operator<< (std::ostream& stream, const Model& model);

   Model();
   Model(size_t sample_size);

   void init();
   
   ~Model();
   
   // Default values;
   const static size_t default_pop_size = 10000;
   const double default_growth_rate = 0.0;
   const double default_mig_rate = 0.0;
   //const double default_growth_rate;
   //const double default_mig_rate;

   // Getters

   double mutation_rate() const { return mutation_rate_; }
   double recombination_rate() const { return recombination_rate_; }
   size_t loci_length() const { return loci_length_; }
   size_t mutation_exact_number() const { return mutation_exact_number_; };

   void set_mutation_rate(double rate) { mutation_rate_ = rate; }
   void set_recombination_rate(double rate) { recombination_rate_ = rate; }
   void set_loci_length(size_t length) { loci_length_ = length; }
   void set_mutation_exact_number(size_t number) { mutation_exact_number_ = number; }

   double growth_rate(size_t pop = 0) const {
     if (current_growth_rates_ == NULL) return default_growth_rate;
     return current_growth_rates_->at(pop);
   }
   
   double population_size(size_t pop = 0) const { 
     if (current_pop_sizes_ == NULL) return default_pop_size;
     return current_pop_sizes_->at(pop);
   }
   

   /**
    * @brief Returns the current migration rate for a given pair of populations.
    *
    * @param source The population from which the migrants come (when looking
    *               backwards in time!)
    * @param sink The population that the migration goes to.
    *
    * @return The current unscaled migration rate.
    */
   double migration_rate(const size_t &source, const size_t &sink) const {
     if (current_mig_rates_ == NULL) return default_mig_rate;
     if (sink == source) return 0.0;
     return current_mig_rates_->at( getMigMatrixIndex(source, sink) );  
   };

   double total_migration_rate(const size_t &sink) const {
     if (current_total_mig_rates_ == NULL) return default_mig_rate;
     return current_total_mig_rates_->at(sink); 
   }; 

   double single_mig_pop(const size_t &source, const size_t &sink) const {
    if (single_mig_probs_list_.at(current_time_idx_) == NULL) return 0.0;
    if (sink == source) return 0.0;
    return single_mig_probs_list_.at(current_time_idx_)->at( getMigMatrixIndex(source, sink) ); 
   }

   bool hasFixedTimeEvent(const double &at_time) const {
    if (single_mig_probs_list_.at(current_time_idx_) == NULL) return false; 
    if (getCurrentTime() != at_time) return false;
    return true;
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
   void set_population_number(const size_t &pop_number) { 
    pop_number_ = pop_number; 
    if (pop_number_<1) throw std::out_of_range("Population number out of range"); 
   }

   void resetTime() { 
     current_pop_sizes_ = pop_sizes_list_.at(0);
     current_growth_rates_ = growth_rates_list_.at(0);
     current_mig_rates_ = mig_rates_list_.at(0);
     current_total_mig_rates_ = total_mig_rates_list_.at(0);
     current_time_idx_ = 0;
   };

   void increaseTime() { 
     if ( current_time_idx_ == change_times_.size() - 1) throw std::out_of_range("Model change times out of range");
     ++current_time_idx_;

     if ( pop_sizes_list_.at(current_time_idx_) != NULL ) 
       current_pop_sizes_ = pop_sizes_list_.at(current_time_idx_);
     if ( growth_rates_list_.at(current_time_idx_) != NULL ) 
       current_growth_rates_ = growth_rates_list_.at(current_time_idx_); 
     if ( mig_rates_list_.at(current_time_idx_) != NULL ) 
       current_mig_rates_ = mig_rates_list_.at(current_time_idx_); 
     if ( total_mig_rates_list_.at(current_time_idx_) != NULL ) 
       current_total_mig_rates_ = total_mig_rates_list_.at(current_time_idx_); 
   };
  
   void print(std::ostream &os) const;


   size_t loci_number() const { return loci_number_; };
   void set_loci_number(size_t loci_number) { loci_number_ = loci_number; }; 
  
   // Add populations size changes
   void addPopulationSizes(double time, const std::vector<double> &pop_sizes, bool relative = false);
   void addPopulationSizes(const double &time, const double &pop_size, bool relative = false); 
   void addPopulationSize(const double &time, const size_t &pop, const double &population_sizes, bool relative = false);

   // Add exponential growth
   void addGrowthRates(const double &time, const std::vector<double> &growth_rates);
   void addGrowthRates(const double &time, const double &growth_rates);
   void addGrowthRate(const double &time, const size_t &population, const double &growth_rates);

   void addSampleSizes(double time, const std::vector<size_t> &samples_sizes);

   // functions to add Migration
   void addMigrationRates(double time, const std::vector<double> &mig_rates);
   void addMigrationRate(double time, size_t source, size_t sink, double mig_rate); 
   void addSymmetricMigration(const double &time, const double &mig_rate); 
   void addSingleMigrationEvent(const double &time, const size_t &source_pop, 
                                const size_t &sink_pop, const double &fraction);

   void finalize(); 

  private:
   Model(const Model&);
   Model& operator=(const Model&);

   size_t addChangeTime(double time);
   
   template <typename T>
   void deleteParList(std::vector<T*> &parList) {
    typename std::vector<T*>::iterator it;
    for (it = parList.begin(); it != parList.end(); ++it) {
      if (*it != NULL) delete *it;
    }
    parList.clear();
   }

   void updateTotalMigRates(const size_t &position);

  void fillVectorList(std::vector<std::vector<double>*> &vector_list, const double &default_value);

  size_t getMigMatrixIndex(const size_t &i, const size_t &j) const {
    assert(i != j);
    return i * (population_number()-1) + j - ( i < j );
  }

   std::vector<size_t> sample_populations_;
   std::vector<double> sample_times_;

   std::vector<double> change_times_;
   std::vector<std::vector<double>*> pop_sizes_list_;
   std::vector<std::vector<double>*> growth_rates_list_;
   std::vector<std::vector<double>*> mig_rates_list_;
   std::vector<std::vector<double>*> total_mig_rates_list_;
   std::vector<std::vector<double>*> single_mig_probs_list_;
   
   size_t current_time_idx_;
   std::vector<double>* current_pop_sizes_;
   std::vector<double>* current_growth_rates_;
   std::vector<double>* current_mig_rates_;
   std::vector<double>* current_total_mig_rates_;

   double mutation_rate_;
   size_t mutation_exact_number_;
   double recombination_rate_;
   size_t pop_number_;

   size_t loci_number_;
   size_t loci_length_;

   size_t exact_window_length_;
   size_t prune_interval_;
};

std::ostream& operator<<(std::ostream& os, const Model& model); 

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &vec) {
  typename std::vector<T>::const_iterator it;
  for (it = vec.begin(); it != vec.end(); ++it) os << *it << " ";
  return os;
}

#endif

/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 *
 * Copyright (C) 2013, 2014 Paul R. Staab, Sha (Joe) Zhu, Dirk Metzler and Gerton Lunter
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

/*!
 * \file model.h
 *
 * \brief This file contains the class Model, which is a simple container object for
 * model parameters.
 *
 */

#ifndef scrm_src_model
#define scrm_src_model
#pragma GCC diagnostic ignored "-Wsign-compare"

#include "macros.h" // Needs to be before cassert

#include <cstddef>
#include <vector>
#include <cfloat>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>

#include "summary_statistics/summary_statistic.h"

class Param;

enum SeqScale { relative, absolute, ms };

class Model
{
  public:
#ifdef UNITTEST
   friend class TestModel;
   friend class TestTimeInterval;
   friend class TestParam;
   friend class TestForest;
   friend class TestPfParam;
#endif
   friend class Forest;
   friend class ForestState;
   friend class Param;
   friend class PfParam;
   friend class CountModel;
   friend std::ostream& operator<< (std::ostream& stream, const Model& model);

   Model();
   Model(size_t sample_size);


   // Default values;
   constexpr static double default_pop_size_ = 10000.0;
   constexpr static double default_growth_rate = 0.0;
   constexpr static double default_mig_rate = 0.0;
   constexpr static double scaling_factor_ = 1.0 / (4 * default_pop_size_);

   // Getters & Setters
   double default_pop_size() const { return Model::default_pop_size_; };


    std::vector<double> change_times() const { return change_times_; }

   /**
    * @brief Returns the scaling factor for times and many parameters
    *
    * @return 1 / ( 4 * default_pop_size);
    */
   double scaling_factor() const { return scaling_factor_; };

   /**
    * @brief Returns the mutation rate per base pair per generation for the
    * currently active sequence position.
    *
    * @return The mutation rate per base per generation
    */
   double mutation_rate() const { return mutation_rates_.at(current_seq_idx_); }

   /**
    * @brief Returns the recombination rate per base pair per generation for the
    * currently active sequence positions.
    *
    * @return The recombination rate per base pair per generation
    */
   double recombination_rate(const size_t idx = -1) const {
     if (idx == -1) return recombination_rates_.at(current_seq_idx_);
     else return recombination_rates_.at(idx);
   }

   /**
    * @brief Returns if the model has recombination.
    *
    * @return true if the model has recombination, false otherwise
    */
   bool has_recombination() const {
     return has_recombination_;
   };


   /**
    * @brief Returns the length of all loci, in base pairs
    *
    * @return length of all loci, in base pairs
    */
   size_t loci_length() const { return loci_length_; }

   /**
    * @brief Getter for the current growth rate of a subpopulation
    *
    * This returns the growth rate for a population for the current time of the
    * model. Use resetTime() and increaseTime() to set the model to the time you
    * want.
    *
    * @param pop The population for which the growth rate is returned.
    *
    * @return The growth rate.
    */
   double growth_rate(size_t pop = 0) const {
     if (current_growth_rates_ == NULL) return default_growth_rate;
     return current_growth_rates_->at(pop);
   }


   /**
    * @brief Getter for the current diploid size of a (sub-)population.
    *
    * This returns the size of a population for the current time of the
    * model. Use resetTime() and increaseTime() to set the model to the time you
    * want. In case of a growing population, you can use the time parameter to
    * specify for which time of the current model stop you want to get the
    * population size.
    *
    * @param pop The population for which the size is returned.
    * @param time The time inside the current model step for which to return the
    * size. This will only make a difference if the model is growing.
    *
    * @return The size of the sub population
    */
   double population_size(const size_t pop = 0, const double time = -1) const {
     return 0.5 / inv_double_pop_size(pop, time);
   }

   double inv_double_pop_size(const size_t pop = 0, const double time = -1) const {
     double pop_size;
     if (current_pop_sizes_ == NULL)
       pop_size = 1/(2*default_pop_size());
     else
       pop_size = current_pop_sizes_->at(pop);   // population size basal to this segment

     if (time >= 0 && growth_rate(pop) != 0.0) {
       assert( time >= getCurrentTime() && time <= getNextTime() );
       pop_size *= std::exp(growth_rate(pop) * (time - getCurrentTime()));
     }

     return pop_size;
   }

   /**
    * @brief Returns the current migration rate for a given pair of populations.
    *
    * The migration rate is returned unscaled, e.g. in migrants per generation.
    * The rate is for the current time of the
    * model. Use resetTime() and increaseTime() to set the model to the time you
    * want.
    *
    * @param source The population from which the migrants come (when looking
    *               backwards in time!)
    * @param sink The population that the migration goes to.
    *
    * @return The current unscaled, backwards migration rate.
    */
   double migration_rate(const size_t source, const size_t sink) const {
     if (sink == source) return 0.0;
     if (current_mig_rates_ == NULL) return default_mig_rate;
     return current_mig_rates_->at( getMigMatrixIndex(source, sink) );
   };

   /**
    * @brief Getter for the current total rate of migrating out of one
    * population (viewed backwards in time).
    *
    *
    * The migration rate is returned unscaled, e.g. in migrants per generation.
    * The rate is for the current time of the
    * model. Use resetTime() and increaseTime() to set the model to the time you
    * want.
    *
    * @param source The population for which the rate is given.
    *
    * @return The total current, unscaled rate of migration out if the population.
    */
   double total_migration_rate(const size_t source) const {
     if (current_total_mig_rates_ == NULL) return default_mig_rate;
     return current_total_mig_rates_->at(source);
   };

   /**
    * @brief Getter for the probability of spontaneous migration at the
    * beginning of the current time interval.
    *
    * Use resetTime() and increaseTime() to set the model to the time interval you
    * want. You can use hasFixedTimeEvent() to check if there is any single migration event.
    *
    * @param source The source population of the migration.
    * @param sink The sink population of the migration.
    *
    * @return The probability/fraction of migration.
    */
   double single_mig_pop(const size_t source, const size_t sink) const {
    if (single_mig_probs_list_.at(current_time_idx_).empty()) return 0.0;
    if (sink == source) return 0.0;
    return single_mig_probs_list_.at(current_time_idx_).at( getMigMatrixIndex(source, sink) );
   }

   void setMutationRate(double rate,
                        const bool &per_locus = false,
                        const bool &scaled = false,
                        const double seq_position = 0);

   void setRecombinationRate(double rate,
                             const bool &per_locus = false,
                             const bool &scaled = false,
                             const double seq_position = 0);

   bool hasFixedTimeEvent(const double at_time) const {
     if (single_mig_probs_list_.at(current_time_idx_).empty()) return false;
     if (getCurrentTime() != at_time) return false;
     return true;
   }

   size_t sample_size() const { return sample_times_.size(); };
   size_t sample_population(size_t sample_id) const { return sample_populations_.at(sample_id); };
   double sample_time(size_t sample_id) const { return sample_times_.at(sample_id); };

   size_t population_number() const { return pop_number_; }

   double getCurrentTime() const { return change_times_.at(current_time_idx_); }
   double getNextTime() const {
    if ( current_time_idx_ + 1 >= change_times_.size() ) return DBL_MAX;
    else return change_times_.at(current_time_idx_ + 1);
   }
   size_t getTimeIdx( double time ) const { // time intervals are half-open [start,end).  First interval is of form [0,end)
      assert (time >= 0);
      size_t idx = 0;
      for (;idx < change_times_.size(); idx++) {
          if (time < change_times_.at(idx)) break;
      }
      return idx - 1;
   }

   double getCurrentSequencePosition() const {
    // Returns position of the next (recombination) rate change
    if ( current_seq_idx_ >= change_position_.size() ) return loci_length();
    return change_position_.at(current_seq_idx_);
   }

   double getNextSequencePosition() const {
    if ( current_seq_idx_ + 1 >= change_position_.size() ) return loci_length();
    else return change_position_.at(current_seq_idx_ + 1);
   }

   double window_length_seq() const { return window_length_seq_; }
   size_t window_length_rec() const { return window_length_rec_; }
   bool has_window_rec() const { return has_window_rec_; }
   bool has_window_seq() const { return has_window_seq_; }
   bool has_approximation() const { return has_appr_; }
   void set_window_length_seq(const double ewl) {
     if (ewl < 0) throw std::invalid_argument("Exact window length can not be negative");
     window_length_seq_ = ewl;
     has_window_seq_ = true;
     has_window_rec_ = false;
     has_appr_ = true;
   }
   void set_window_length_rec(const size_t ewl) {
     window_length_rec_ = ewl;
     has_window_seq_ = false;
     has_window_rec_ = true;
     has_appr_ = true;
   }
   void disable_approximation() {
     has_appr_ = false;
     has_window_rec_ = false;
     has_window_seq_ = false;
   }

   void set_population_number(const size_t pop_number) {
    pop_number_ = pop_number;
    if (pop_number_<1) throw std::out_of_range("Population number out of range");
   }

   void resetTime( ) {
     if (pop_sizes_list_[0].empty()) current_pop_sizes_ = NULL;
     else current_pop_sizes_ = &(pop_sizes_list_[0]);
     if (growth_rates_list_[0].empty()) current_growth_rates_ = NULL;
     else current_growth_rates_ = &(growth_rates_list_[0]);
     if (mig_rates_list_[0].empty()) current_mig_rates_ = NULL;
     else current_mig_rates_ = &(mig_rates_list_[0]);
     if (total_mig_rates_list_[0].empty()) current_total_mig_rates_ = NULL;
     else current_total_mig_rates_ = &(total_mig_rates_list_[0]);
     current_time_idx_ = 0;
   };

   void resetTime( double current_time ) {
     current_time_idx_ = 0;
     while (getNextTime() <= current_time ) {
        if ( current_time_idx_ == change_times_.size() - 1) throw std::out_of_range("Model change times out of range");
        ++current_time_idx_;
     }
     // is the check for non-NULL-ness redundant?  Change into an assert?
     if ( !pop_sizes_list_.at(current_time_idx_).empty() )
       current_pop_sizes_ = &(pop_sizes_list_.at(current_time_idx_));
     if ( !growth_rates_list_.at(current_time_idx_).empty() )
       current_growth_rates_ = &(growth_rates_list_.at(current_time_idx_));
     if ( !mig_rates_list_.at(current_time_idx_).empty() )
       current_mig_rates_ = &(mig_rates_list_.at(current_time_idx_));
     if ( !total_mig_rates_list_.at(current_time_idx_).empty() )
       current_total_mig_rates_ = &(total_mig_rates_list_.at(current_time_idx_));
   }

   void resetSequencePosition() {
     current_seq_idx_ = 0;
   }

   void increaseTime() {
     if ( current_time_idx_ == change_times_.size() ) throw std::out_of_range("Model change times out of range");
     ++current_time_idx_;

     if ( ! pop_sizes_list_.at(current_time_idx_).empty() )
       current_pop_sizes_ = &(pop_sizes_list_.at(current_time_idx_));
     if ( ! growth_rates_list_.at(current_time_idx_).empty() )
       current_growth_rates_ = &(growth_rates_list_.at(current_time_idx_));
     if ( ! mig_rates_list_.at(current_time_idx_).empty() )
       current_mig_rates_ = &(mig_rates_list_.at(current_time_idx_));
     if ( ! total_mig_rates_list_.at(current_time_idx_).empty() )
       current_total_mig_rates_ = &(total_mig_rates_list_.at(current_time_idx_));
//=======  JOE: I think the following can be replaced
     //// is the check for non-NULL-ness redundant?  Change into an assert?
     //if ( pop_sizes_list_.at(current_time_idx_) != NULL )
       //current_pop_sizes_ = pop_sizes_list_.at(current_time_idx_);
     //if ( growth_rates_list_.at(current_time_idx_) != NULL )
       //current_growth_rates_ = growth_rates_list_.at(current_time_idx_);
     //if ( mig_rates_list_.at(current_time_idx_) != NULL )
       //current_mig_rates_ = mig_rates_list_.at(current_time_idx_);
     //if ( total_mig_rates_list_.at(current_time_idx_) != NULL )
       //current_total_mig_rates_ = total_mig_rates_list_.at(current_time_idx_);
//>>>>>>> squashing_JZ
   };

   size_t getNumEpochs() const { return change_times_.size(); }

   void increaseSequencePosition() {
    ++current_seq_idx_;
   }

   size_t countChangeTimes() const { return change_times_.size(); }
   size_t countChangePositions() const { return change_position_.size(); }
   void decreaseSequencePosition() {
    --current_seq_idx_;
   }

   void print(std::ostream &os) const;


   size_t loci_number() const { return loci_number_; };
   void set_loci_number(size_t loci_number) { loci_number_ = loci_number; };

   // Add populations size changes
   void addPopulationSizes(double time, const std::vector<double> &pop_sizes,
                           const bool &time_scaled = false, const bool &relative = false);

   void addPopulationSizes(const double time, const double pop_size,
                           const bool &time_scaled = false, const bool &relative = false);

   void addPopulationSize(const double time, const size_t pop, double population_sizes,
                          const bool &time_scaled = false, const bool &relative = false);

   // Add exponential growth
   void addGrowthRates(const double time, const std::vector<double> &growth_rates,
                       const bool &time_scaled = false,
                       const bool &rate_scaled = false);

   void addGrowthRates(const double time, const double growth_rates,
                       const bool &time_scaled = false,
                       const bool &rate_scaled = false);

   void addGrowthRate(const double time, const size_t population,
                      double growth_rates, const bool &time_scaled = false,
                      const bool &rate_scaled = false);

   void addSampleSizes(double time, const std::vector<size_t> &samples_sizes,
                         const bool &scaled = false);

   // functions to add Migration
   void addMigrationRates(double time, const std::vector<double> &mig_rates,
                          const bool &time_scaled = false, const bool &rate_scaled = false);

   void addMigrationRate(double time, size_t source, size_t sink, double mig_rate,
                         const bool &time_scaled = false, const bool &rate_scaled = false);

   void addSymmetricMigration(const double time, const double mig_rate,
                              const bool &time_scaled = false, const bool &rate_scaled = false);

   void addSingleMigrationEvent(const double time, const size_t source_pop,
                                const size_t sink_pop, const double fraction,
                                const bool &time_scaled = false);

   void finalize();

   void check();
   //void reset();

   size_t countSummaryStatistics() const {
     return summary_statistics_.size();
   }

   SummaryStatistic* getSummaryStatistic(const size_t i) const {
     return summary_statistics_.at(i).get();
   }

   void addSummaryStatistic(std::shared_ptr<SummaryStatistic> sum_stat) {
     summary_statistics_.push_back(sum_stat);
   }

  void addPopulation();

  SeqScale getSequenceScaling() const { return seq_scale_; }
  void setSequenceScaling(SeqScale seq_scale) { seq_scale_ = seq_scale; };

   void setLocusLength(const size_t length) {
    // Rescale the rates that are per base pair
    for (size_t i = 0; i < change_position_.size(); ++i) {
      mutation_rates_.at(i) *= (double)loci_length() / length;
      recombination_rates_.at(i) *= (double)(loci_length()-1) / (length-1);
    }
    loci_length_ = length;
   }

   //biased sampling
    bool biased_sampling = false;                         // this is turned on if br or bh is set
    std::vector<double> bias_heights_;                    // this vector defines the outer time boundaries of the sections with different rec rate
    std::vector<double> bias_strengths_;                  // this vector defines how the below ratios should be scaled relative to each other
    std::vector<double> bias_ratios_;                     // this vector holds ratios r, to impose rec rate of rho*r_i in time section i
    // NB: bias_ratios_ must be scaled for each demographic model using model_summary, otherwise the recombination rate estimate will diverge
    std::vector<double> application_delays;                 // should init to (getNumEpochs(), 0)

    const std::vector<double>& bias_heights() const {return bias_heights_;}
    const std::vector<double>& bias_ratios() const {return bias_ratios_;}
    const std::vector<double>& bias_strengths() const {return bias_strengths_;}

    void clearBiasHeights()                  { bias_heights_.clear(); }
    void clearBiasRatios()                   { bias_ratios_.clear(); }
    void clearBiasStrengths()                { bias_strengths_.clear(); }
    void addToBiasHeights(double height)     { biased_sampling = true;
                                               bias_heights_.push_back( height ); } //height is scaled in generations
    void addToBiasRatios(double ratio)       { biased_sampling = true;
                                               bias_ratios_.push_back( ratio ); }
    void addToBiasStrengths(double strength) { biased_sampling = true;
                                               bias_strengths_.push_back( strength ); }

    void lags_to_application_delays(std::vector<double> lags) {
        application_delays.resize(0);
        for( size_t i=0; i<lags.size(); i++){
            application_delays.push_back( 0.5*lags.at(i) );
        }
    }

    const std::vector<std::vector<double> >& single_mig_probs_list() const { return single_mig_probs_list_; }

  private:

   std::vector<double> change_times_;

   double change_position(size_t idx) const {
    return this->change_position_.at(idx);
   }

   size_t get_position_index() const {
    return this->current_seq_idx_;
   }

   size_t addChangeTime(double time, const bool &scaled = false);
   size_t addChangePosition(const double position);

   template <typename T>
   void deleteParList(std::vector<T*> &parList) {
    typename std::vector<T*>::iterator it;
    for (it = parList.begin(); it != parList.end(); ++it) {
      if (*it != NULL) delete *it;
    }
    parList.clear();
   }

   void updateTotalMigRates(const size_t position);
   bool has_migration() { return has_migration_; };

  void fillVectorList(std::vector<std::vector<double> > &vector_list, const double default_value);
  void calcPopSizes();
  void checkPopulation(const size_t pop) {
    if (pop >= this->population_number())
      throw std::invalid_argument("Invalid population");
  }

  friend void swap(Model& first, Model& second);

  size_t getMigMatrixIndex(const size_t i, const size_t j) const {
    assert(i != j);
    return i * (population_number()-1) + j - ( i < j );
  }

  void addPopToMatrixList(std::vector<std::vector<double> > &vector_list,
                          size_t new_pop,
                          double default_value = nan("value to replace"));
  void addPopToVectorList(std::vector<std::vector<double> > &vector_list);

   // Stores information about samples. Each index represents a sample.
   std::vector<size_t> sample_populations_;
   std::vector<double> sample_times_;

   // Stores the time and sequences positions where the model changes.
   std::vector<double> change_position_;

   // These pointer vectors hold the actual model parameters that can change in
   // time. Each index represents one period in time within which the model
   // parameters are constant. NULL means that the parameters do not change.
   std::vector<std::vector<double> > growth_rates_list_;
   std::vector<std::vector<double> > mig_rates_list_;
   std::vector<std::vector<double> > total_mig_rates_list_;
   std::vector<std::vector<double> > single_mig_probs_list_;

   // Population sizes are saved as 1/(2N), where N is the actual population
   // size (do to fast multiplication rather than slow division in the
   // algorithm)
   std::vector<std::vector<double> > pop_sizes_list_;

   // These vectors contain the model parameters that may change along the sequence.
   // Again, each index represents an sequence segment within with the model
   // does not change.
   std::vector<double> recombination_rates_;       /*!< Unit: Recombinations per base per generation */
   std::vector<double> mutation_rates_;           /*!< Unit: Mutations per base per generation */

   // The index of the time and sequence segment currently active.
   size_t current_time_idx_;
   size_t current_seq_idx_;

   // Direct pointers to the currently active model parameters.
   std::vector<double>* current_pop_sizes_;
   std::vector<double>* current_growth_rates_;
   std::vector<double>* current_mig_rates_;
   std::vector<double>* current_total_mig_rates_;

   size_t pop_number_;

   size_t loci_number_;
   size_t loci_length_;

   double window_length_seq_;
   size_t window_length_rec_;
   bool has_window_seq_;
   bool has_window_rec_;
   bool has_appr_;

   bool has_migration_;
   bool has_recombination_;

   SeqScale seq_scale_;

   std::vector<std::shared_ptr<SummaryStatistic> > summary_statistics_;


};



std::ostream& operator<<(std::ostream& os, Model& model);

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &vec) {
  typename std::vector<T>::const_iterator it;
  for (it = vec.begin(); it != vec.end(); ++it) os << *it << " ";
  return os;
}



#endif

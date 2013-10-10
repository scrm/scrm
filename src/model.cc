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

#include "model.h"


Model::Model() { 
  this->init();
}

Model::Model(size_t sample_size) {
  this->init();

  this->addSampleSizes(0.0, std::vector<size_t>(1, sample_size));
  this->resetTime();
}

Model::~Model() { 
  deleteParList(pop_sizes_list_);
  deleteParList(growth_rates_list_);
  deleteParList(mig_rates_list_);
  deleteParList(total_mig_rates_list_);
  deleteParList(single_mig_probs_list_);
}

// Copy constructor
Model::Model(const Model& model) {
  //normal members
  default_pop_size = model.default_pop_size;
  default_loci_length = model.default_loci_length;
  default_growth_rate = model.default_growth_rate;
  default_mig_rate = model.default_mig_rate;
  
  mutation_rate_ = model.mutation_rate_;
  mutation_rate_per_locus_ = model.mutation_rate_per_locus_;
  mutation_exact_number_ = model.mutation_exact_number_;
  rec_rate_ = model.rec_rate_;
  pop_number_ = model.pop_number_;
  loci_number_ = model.loci_number_;
  loci_length_ = model.loci_length_;
  exact_window_length_ = model.exact_window_length_;
  prune_interval_ = model.prune_interval_;
  has_migration_ = model.has_migration_;

  // Vector members
  sample_times_ = model.sample_times_;
  sample_populations_ = model.sample_populations_;
  change_times_ = model.change_times_;

  // Vector lists
  pop_sizes_list_ = copyVectorList(model.pop_sizes_list_);
  growth_rates_list_ = copyVectorList(model.growth_rates_list_);
  mig_rates_list_ = copyVectorList(model.mig_rates_list_);
  total_mig_rates_list_ = copyVectorList(model.total_mig_rates_list_);
  single_mig_probs_list_ = copyVectorList(single_mig_probs_list_);

  // Pointers
  current_time_idx_ = model.current_time_idx_; 
  current_pop_sizes_ = model.pop_sizes_list_.at(current_time_idx_);
  current_growth_rates_ = model.growth_rates_list_.at(current_time_idx_); 
  current_mig_rates_ = model.mig_rates_list_.at(current_time_idx_);
  current_total_mig_rates_ = model.total_mig_rates_list_.at(current_time_idx_);
}

void Model::init() {
  default_pop_size = 10000;
  default_loci_length = 100000;
  default_growth_rate = 0.0;
  default_mig_rate = 0.0;

  sample_times_ = std::vector<double>();
  sample_populations_ = std::vector<size_t>();

  change_times_ = std::vector<double>();
  pop_sizes_list_ = std::vector<std::vector<double>*>();
  growth_rates_list_ = std::vector<std::vector<double>*>();
  mig_rates_list_ = std::vector<std::vector<double>*>();
  total_mig_rates_list_ = std::vector<std::vector<double>*>(); 
  single_mig_probs_list_ = std::vector<std::vector<double>*>();

  has_migration_ = false;

  this->addChangeTime(0.0);

  this->set_population_number(1);

  this->set_loci_number(1);
  this->loci_length_ = this->default_loci_length;

  this->set_mutation_rate(0.0);
  this->set_recombination_rate(0.0, default_loci_length);

  this->set_exact_window_length(-1);
  this->set_prune_interval(10);

  this->mutation_exact_number_ = -1;

  this->resetTime();
}


void Model::reset() {
  deleteParList(pop_sizes_list_);
  deleteParList(growth_rates_list_);
  deleteParList(mig_rates_list_);
  deleteParList(total_mig_rates_list_);
  deleteParList(single_mig_probs_list_);

  init();
}

/** 
 * Function to add a new change time to the model.
 *
 * It preserves the relation between the times and the *param*_list_ containers.
 * If the same time is added multiple times, it is just added once to the model,
 * but this should not make a difference when using this function.
 *
 * @param time The time that is added
 * @param scaled set to TRUE if the time is in units of 4N0 generations, and
 * FALSE if it is in units of generations. 
 *
 * @returns The position the time has now in the vector
 */
size_t Model::addChangeTime(double time, const bool &scaled) {
  if (scaled) time *= 4 * default_pop_size;

  size_t position = 0;
  if ( change_times_.size() == 0 ) {
    change_times_ = std::vector<double>(1, time);
    pop_sizes_list_.push_back(NULL);
    growth_rates_list_.push_back(NULL);
    mig_rates_list_.push_back(NULL);
    total_mig_rates_list_.push_back(NULL);
    single_mig_probs_list_.push_back(NULL);
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
  pop_sizes_list_.insert(pop_sizes_list_.begin() + position, NULL);
  growth_rates_list_.insert(growth_rates_list_.begin() + position, NULL);
  mig_rates_list_.insert(mig_rates_list_.begin() + position, NULL);
  total_mig_rates_list_.insert(total_mig_rates_list_.begin() + position, NULL);
  single_mig_probs_list_.insert(single_mig_probs_list_.begin() + position, NULL);
  return position;
}


void Model::addSampleSizes(double time, const std::vector<size_t> &samples_sizes) {
  for (size_t pop = 0; pop < samples_sizes.size(); ++pop) {
    for (size_t i = 0; i < samples_sizes.at(pop); ++i) {
      sample_populations_.push_back(pop);
      sample_times_.push_back(time);
    }
  }
}


/**
 * @brief This changes the size of all populations to the given values at a
 * specific point in time. 
 *
 * The sizes apply for with point on backwards in time.    
 *
 * @param time The time at which the population changes their sizes.
 * @param pop_sizes A vector of new population sizes. Can either be given as
 *    fraction of N0 or as an absolute value. See relative.
 * @param time_scaled Set to true if the time is given in units of 4*N0
 *    generations, or to false if the time is given in units of generations.
 * @param relative set to TRUE, if the population sizes are given relative to
 *    N0, or to FALSE if they are absolute values.
 */
void Model::addPopulationSizes(double time, const std::vector<double> &pop_sizes, 
                               const bool &time_scaled, const bool &relative) {

  if ( pop_sizes.size() != population_number() ) 
    throw std::logic_error("Population size values do not meet the number of populations");
  auto pop_sizes_heap = new std::vector<double>(pop_sizes);
  if (relative) {
    for (auto it = pop_sizes_heap->begin(); it != pop_sizes_heap->end(); ++it) {
      if (std::isnan(*it)) continue;
      else *it *= this->default_pop_size; 
    }
  }
  size_t position = addChangeTime(time, time_scaled);
  pop_sizes_list_[position] = pop_sizes_heap;  
}


/**
 * @brief This changes the size of all populations to a given value at a
 * specific point in time. 
 *
 * The sizes apply for with point on backwards in time.    
 *
 * @param time The time at which the population changes their sizes.
 * @param pop_sizes The size to which we set all populations. Can either be given as
 *    fraction of N0 or as an absolute value. See relative.
 * @param time_scaled Set to true if the time is given in units of 4*N0
 *    generations, or to false if the time is given in units of generations.
 * @param relative set to TRUE, if the population sizes are given relative to
 *    N0, or to FALSE if they are absolute values.
 */
void Model::addPopulationSizes(const double &time, const double &pop_size, 
                               const bool &time_scaled, const bool &relative) {
  addPopulationSizes(time, std::vector<double>(population_number(), pop_size), time_scaled, relative);
}
  

/**
 * @brief This changes the size of a single populations to a given value at a
 * specific point in time. 
 *
 * The sizes apply for with point on backwards in time.    
 * Requires Model.finalization() to be called after the model is set up.
 *
 * @param time The time at which the population change its size.
 * @param pop The population which will change its size.
 * @param population_size The size to which we set the population. Can either be given as
 *    fraction of N0 or as an absolute value. See relative.
 * @param time_scaled Set to true if the time is given in units of 4*N0
 *    generations, or to false if the time is given in units of generations.
 * @param relative set to TRUE, if the population sizes are given relative to
 *    N0, or to FALSE if they are absolute values.
 */
void Model::addPopulationSize(const double &time, const size_t &pop, double population_size,
                              const bool &time_scaled, const bool &relative) {
  size_t position = addChangeTime(time, time_scaled);
  if (relative) population_size *= default_pop_size;

  if (pop_sizes_list_.at(position) == NULL) addPopulationSizes(time, nan("value to replace"), time_scaled);
  pop_sizes_list_.at(position)->at(pop) = population_size;
}


/**
 * @brief Set the population size growth rates at a certain time point.
 *
 * The population growth or shrinks exponentially from that time point on
 * backwards in time. 
 * Requires Model.finalization() to be called after the model is set up.
 *
 * @param time The time at which to set the growth rates
 * @param growth_rates A vector of growth rates for all populations
 * @param time_scaled Set to true if the time is given in units of 4*N0
 *    generations, or to false if the time is given in units of generations.
 */
void Model::addGrowthRates(const double &time, const std::vector<double> &growth_rates,
                           const bool &time_scaled) {
  if ( growth_rates.size() != population_number() ) 
    throw std::logic_error("Growth rates values do not meet the number of populations");
  std::vector<double>* growth_rates_heap = new std::vector<double>(growth_rates);
  size_t position = addChangeTime(time, time_scaled);
  growth_rates_list_[position] = growth_rates_heap; 

}


/**
 * @brief Set the population size growth rates at a certain time point.
 *
 * The population growth or shrinks exponentially from that time point on
 * backwards in time. 
 * Requires Model.finalization() to be called after the model is set up.
 *
 * @param time The time at which to set the growth rates
 * @param growth_rates The growth rate for all populations
 * @param time_scaled Set to true if the time is given in units of 4*N0
 *    generations, or to false if the time is given in units of generations.
 */
void Model::addGrowthRates(const double &time, const double &growth_rate,
                           const bool &time_scaled) {
  addGrowthRates(time, std::vector<double>(population_number(), growth_rate), time_scaled);
}


/**
 * @brief Set the population size growth rates of a population at a certain time point.
 *
 * The population growth or shrinks exponentially from that time point on
 * backwards in time. 
 * Requires Model.finalization() to be called after the model is set up.
 *
 * @param time The time at which to set the growth rates
 * @param population The population to which the growth rate applies.
 * @param growth_rates The growth rate for the populations
 * @param time_scaled Set to true if the time is given in units of 4*N0
 *    generations, or to false if the time is given in units of generations.
 */
void Model::addGrowthRate(const double &time, const size_t &population, 
                          const double &growth_rate, const bool &time_scaled) {
  size_t position = addChangeTime(time, time_scaled);
  if (growth_rates_list_.at(position) == NULL) addGrowthRates(time, nan("number to replace"), time_scaled); 
  growth_rates_list_.at(position)->at(population) = growth_rate;
}


/**
 * @brief Sets a migration rate form a specific population to another starting from a
 * certain time point (going backwards in time); 
 *
 * This requires model finalization, e.g. call model.finalize() after you set up
 * the model completely.
 *
 * @param time The time at which the migration is set to the given value.
 *        It applies backwards in time until it is changed again.
 * @param source The population from which the individuals migrate from when
 *        looking backwards in time. Is the sink population when looking forward.
 * @param sink The population to which the individuals migrate to (also
 *        when looking backwards in time)
 * @param mig_rate The backwards scaled migration rate M_ij = 4N0 * m_ij, 
 *        where m_ij is the fraction for population i = source that migrates 
 *        to population j = sink (again, when looking backwards in time).
 * @param scaled_time Set to true if the time is given in units of 4*N0
 *    generations, or to false if the time is given in units of generations.
 * @param scaled_rate Set to true if the rate is given as M = 4*N0*m and to
 *  false if it is given as m.
 *
 */
void Model::addMigrationRate(double time, size_t source, size_t sink, double mig_rate,
                             const bool &scaled_time, const bool &scaled_rates) {
  size_t position = addChangeTime(time, scaled_time);
  if (scaled_rates) mig_rate /= 4 * default_pop_size;
  if (mig_rates_list_.at(position) == NULL) addSymmetricMigration(time, nan("value to replace"), scaled_time); 
  mig_rates_list_.at(position)->at(getMigMatrixIndex(source, sink)) = mig_rate;  
}


/** 
 * @brief Sets the migration matrix to the given values for the time following at
 * certain time point (backwards in time). 
 *
 * This requires model finalization, e.g. call model.finalize() after you set up
 * the model completely.
 *
 * @param time The time at which the migration is set to the given values.
 *        The values apply backwards in time until they are changed again.
 * @param mig_rates The (backwards) scaled migration matrix, given as concatenation of
 *        row vectors (M_00, M_01, ..., M_0n, M_10, ..., M_n1, ..., M_nn), where
 *        M_ij = 4N0 * m_ij and m_ij is the faction of population i that
 *        migrates to population j (viewed backwards in time; forwards the
 *        migration is from population j to i). The diagonal elements of the
 *        matrix are ignored and can be set to "x" for better readability. 
 * @param scaled_time Set to true if the time is given in units of 4*N0
 *    generations, or to false if the time is given in units of generations.
 * @param scaled_rate Set to true if the rate is given as M = 4*N0*m and to
 *  false if it is given as m.
 */
void Model::addMigrationRates(double time, const std::vector<double> &mig_rates,
                              const bool &scaled_time, const bool &scaled_rates) {
  double popnr = population_number();
  double scaling = 1;
  if (scaled_rates) scaling = 1 / ( 4 * default_pop_size );
  if ( mig_rates.size() != population_number()*population_number() ) 
    throw std::logic_error("Migration rates values do not meet the number of populations");
  std::vector<double>* mig_rates_heap = new std::vector<double>();
  mig_rates_heap->reserve(popnr*popnr-popnr);
  for (size_t i = 0; i < popnr; ++i) {
    for (size_t j = 0; j < popnr; ++j) {
      if (i == j) continue;
      mig_rates_heap->push_back(mig_rates.at(i*popnr+j) * scaling); 
    }
  }

  size_t position = addChangeTime(time, scaled_time);
  mig_rates_list_[position] = mig_rates_heap; 
}



/**
 * @brief Adds symmetric migration to a model.
 *
 * Sets migration between all population to the given value starting at a
 * certain time point going backwards in time. Unlike 'ms', this uses the actual
 * value provided and does not divide it by #population-1.
 *
 * This requires model finalization, e.g. call model.finalize() after you set up
 * the model completely.
 *
 * @param time The time at which the migration is set to the given value.
 *        It applies backwards in time until it is changed again.
 * @param mig_rate The scaled migration rate M_ij = 4N0 * m_ij that is used
 *        between all populations i and j. m_ij is the fraction of population i 
 *        that migrates to population j.
 * @param time_scaled Set to true if the time is given in units of 4*N0
 *    generations, or to false if the time is given in units of generations.
 * @param rate_scaled Set to true if the rate is given as M = 4*N0*m and to
 *  false if it is given as m.
 */
void Model::addSymmetricMigration(const double &time, const double &mig_rate, 
                                  const bool &time_scaled, const bool &rate_scaled) {
  std::vector<double> mig_rates = std::vector<double>(population_number()*population_number(), mig_rate);
  this->addMigrationRates(time, mig_rates, time_scaled, rate_scaled);
}


void Model::addSingleMigrationEvent(const double &time, const size_t &source_pop, 
                                    const size_t &sink_pop, const double &fraction,
                                    const bool &time_scaled) {
  
  size_t position = addChangeTime(time, time_scaled);
  size_t popnr = population_number();
  
  if ( time < 0.0 ) throw std::invalid_argument("Single migration event: Negative time");
  if ( source_pop >= population_number() ) throw std::invalid_argument("Single migration event: Unknown population");
  if ( sink_pop >= population_number() ) throw std::invalid_argument("Single migration event: Unknown population");
  if ( fraction < 0.0 || fraction > 1.0 ) throw std::invalid_argument("Single migration event: Fraction out of range");

  if ( single_mig_probs_list_.at(position) == NULL ) {
    std::vector<double> *mig_probs = new std::vector<double>(popnr*popnr-popnr, 0.0);
    single_mig_probs_list_.at(position) = mig_probs;
  }

  single_mig_probs_list_.at(position)->at(getMigMatrixIndex(source_pop, sink_pop)) = fraction; 
} 


std::ostream& operator<<(std::ostream& os, const Model& model) {
  os << "---- Model: ------------------------" << std::endl;
  os << "Mutation rate: " << model.mutation_rate() << std::endl;  
  os << "Recombination rate: " << model.recombination_rate() << std::endl;  
  os << "Sample size: " << model.sample_size() << std::endl;  
  
  for (size_t idx = 0; idx < model.change_times_.size(); ++idx) { 
    os << std::endl << "At time " << model.change_times_.at(idx) << ":" << std::endl;  
    if (model.pop_sizes_list_.at(idx) != NULL) {
      os << " Population sizes: " << *(model.pop_sizes_list_.at(idx)) << std::endl;
    }
    if (model.growth_rates_list_.at(idx) != NULL) {
      os << " Growth Rates: " << *(model.growth_rates_list_.at(idx)) << std::endl;
    }
    if (model.mig_rates_list_.at(idx) != NULL) {
      os << " Mig Rates: " << *(model.mig_rates_list_.at(idx)) << std::endl;
    }
    if (model.single_mig_probs_list_.at(idx) != NULL) {
      os << " Single Mig Rates: " << *(model.single_mig_probs_list_.at(idx)) << std::endl;
    }
  }
  os << "------------------------------------" << std::endl;
  return(os);
}


void Model::updateTotalMigRates(const size_t &position) {
  std::vector<double>* mig_rates;
  if ( total_mig_rates_list_.at(position) == NULL ) {
    mig_rates = new std::vector<double>(population_number(), 0.0);
    total_mig_rates_list_.at(position) = mig_rates;
  }
  else {
    mig_rates = total_mig_rates_list_.at(position); 
  }

  for (size_t i = 0; i < population_number(); ++i) {
    for (size_t j = 0; j < population_number(); ++j) {
      if (i == j) continue;
      mig_rates->at(i) += mig_rates_list_.at(position)->at( getMigMatrixIndex(i,j) );
    }
    if (mig_rates->at(i) > 0) has_migration_ = true;
  }
}; 


void Model::finalize() {
  fillVectorList(mig_rates_list_, default_mig_rate);
  fillVectorList(pop_sizes_list_, default_pop_size);
  fillVectorList(growth_rates_list_, default_growth_rate);

  for (size_t j = 0; j < mig_rates_list_.size(); ++j) {
    if (mig_rates_list_.at(j) == NULL) continue;
    updateTotalMigRates(j);
  } 
  
  check();
}


void Model::check() {
  // Sufficient sample size?
  if (sample_size() < 2) throw std::invalid_argument("Sample size needs be to at least 2");

  // Structure without migration?
  if (population_number() > 1 && !has_migration())
    throw std::invalid_argument("Model has multiple population but no migration. Coalescence impossible"); 
}


void Model::fillVectorList(std::vector<std::vector<double>*> &vector_list, const double &default_value) {
  std::vector<double>* last = NULL; 
  std::vector<double>* current = NULL; 
  for (size_t j = 0; j < vector_list.size(); ++j) {
    current = vector_list.at(j);
    if (current == NULL) continue;

    for (size_t i = 0; i < current->size(); ++i) {
      if ( !std::isnan(current->at(i)) ) continue;

      if (last == NULL) (current)->at(i) = default_value;
      else current->at(i) = last->at(i); 
    }
    last = current;
  }
}


std::vector<std::vector<double>*> Model::copyVectorList(const std::vector<std::vector<double>*> &source) {
  auto copy = std::vector<std::vector<double>*>();

  for (size_t j = 0; j < source.size(); ++j) {
    if (source.at(j) == NULL) { 
      copy.push_back(NULL);
      continue;
    }
    copy.push_back(new std::vector<double>(*source.at(j)));
  }

  return copy;
}


void swap(Model& first, Model& second) {
  using std::swap;
  swap(first.default_pop_size, second.default_pop_size);
  swap(first.default_loci_length, second.default_loci_length);
  swap(first.default_growth_rate, second.default_growth_rate);
  swap(first.default_mig_rate, second.default_mig_rate);

  swap(first.mutation_rate_, second.mutation_rate_);
  swap(first.mutation_rate_per_locus_, second.mutation_rate_per_locus_);
  swap(first.mutation_exact_number_, second.mutation_exact_number_);
  swap(first.rec_rate_, second.rec_rate_);
  swap(first.pop_number_, second.pop_number_);
  swap(first.loci_number_, second.loci_number_);
  swap(first.loci_length_, second.loci_length_);
  swap(first.exact_window_length_, second.exact_window_length_);
  swap(first.prune_interval_, second.prune_interval_);
  swap(first.has_migration_, second.has_migration_);

  // Vector members
  swap(first.sample_times_, second.sample_times_);
  swap(first.sample_populations_, second.sample_populations_);
  swap(first.change_times_, second.change_times_);

  // Vector lists
  swap(first.pop_sizes_list_, second.pop_sizes_list_);
  swap(first.growth_rates_list_, second.growth_rates_list_);
  swap(first.mig_rates_list_, second.mig_rates_list_);
  swap(first.total_mig_rates_list_, second.total_mig_rates_list_);
  swap(first.single_mig_probs_list_, second.single_mig_probs_list_);

  // Pointers
  swap(first.current_time_idx_, second.current_time_idx_); 
  swap(first.current_pop_sizes_, second.current_pop_sizes_);
  swap(first.current_growth_rates_, second.current_growth_rates_); 
  swap(first.current_mig_rates_, second.current_mig_rates_);
  swap(first.current_total_mig_rates_, second.current_total_mig_rates_);
}

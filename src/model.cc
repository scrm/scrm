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

#include "model.h"


Model::Model() :
  has_migration_(false),
  has_recombination_(false) {

  this->set_loci_number(1);
  this->setLocusLength(1);
  this->addChangeTime(0.0);
  this->addChangePosition(0.0);

  this->set_population_number(1);

  this->setMutationRate(0.0);
  this->setRecombinationRate(0.0);

  this->window_length_seq_ = 0;
  this->set_window_length_rec(500);

  this->setSequenceScaling(ms);

  this->bias_heights_ = std::vector<double> ({0.0, DBL_MAX});
  this->bias_strengths_ = std::vector <double> ({1.0});

  this->resetTime();
  this->resetSequencePosition();
}


Model::Model(size_t sample_size) :
  has_migration_(false),
  has_recombination_(false) {

  this->set_loci_number(1);
  this->setLocusLength(1);
  this->addChangeTime(0.0);
  this->addChangePosition(0.0);

  this->set_population_number(1);

  this->setMutationRate(0.0);
  this->setRecombinationRate(0.0);

  this->window_length_seq_ = 0;
  this->set_window_length_rec(500);

  this->setSequenceScaling(ms);

  this->bias_heights_ = std::vector<double> ({0.0, DBL_MAX});
  this->bias_strengths_ = std::vector <double> ({1.0});

  this->addSampleSizes(0.0, std::vector<size_t>(1, sample_size));
  this->setLocusLength(1000);
  this->resetTime();
  this->resetSequencePosition();
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
  if (scaled) time *= 4 * default_pop_size();

  size_t position = 0;
  if ( change_times_.size() == 0 ) {
    change_times_ = std::vector<double>(1, time);
    pop_sizes_list_.push_back(std::vector<double>());
    growth_rates_list_.push_back(std::vector<double>());
    mig_rates_list_.push_back(std::vector<double>());
    total_mig_rates_list_.push_back(std::vector<double>());
    single_mig_probs_list_.push_back(std::vector<double>());
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
  pop_sizes_list_.insert(pop_sizes_list_.begin() + position, std::vector<double>());
  growth_rates_list_.insert(growth_rates_list_.begin() + position, std::vector<double>());
  mig_rates_list_.insert(mig_rates_list_.begin() + position, std::vector<double>());
  total_mig_rates_list_.insert(total_mig_rates_list_.begin() + position, std::vector<double>());
  single_mig_probs_list_.insert(single_mig_probs_list_.begin() + position, std::vector<double>());
  return position;
}


/**
 * Adds a new change position to the model.
 *
 * Change position are sequence positions where mutation or recombination rates
 * change. This creates a new position, but does not add the new rates.
 *
 * @param position The sequence position add which a change is added
 *
 * @returns The index of the new rates in the recombination_rates_ and
 * mutation_rates vectors.
 */
size_t Model::addChangePosition(const double position) {
  if (position < 0 || position > loci_length()) {
    std::stringstream ss;
    ss << "Rate change position " << position << " is outside of the simulated sequence.";
    throw std::invalid_argument(ss.str());
  }

  if ( change_position_.size() == 0 ) {
    change_position_ = std::vector<double>(1, position);
    recombination_rates_.push_back(-1);
    mutation_rates_.push_back(-1);
    return 0;
  }

  size_t idx = 0;
  std::vector<double>::iterator ti;
  for (ti = change_position_.begin(); ti != change_position_.end(); ++ti) {
    if ( *ti == position ) return idx;
    if ( *ti > position ) break;
    ++idx;
  }

  change_position_.insert(ti, position);

  // Add Null at the right position in all parameter vectors
  recombination_rates_.insert(recombination_rates_.begin() + idx, -1);
  mutation_rates_.insert(mutation_rates_.begin() + idx, -1);

  return idx;
}


void Model::addSampleSizes(double time, const std::vector<size_t> &samples_sizes, const bool &scaled) {
  if (scaled) time *= 4 * default_pop_size();

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

  size_t position = addChangeTime(time, time_scaled);

  pop_sizes_list_[position].clear();
  for (double pop_size : pop_sizes) {
    if (!std::isnan(pop_size)) {
      // Scale to absolute values if necessary
      if (relative) { pop_size *= this->default_pop_size(); }

      // Save inverse double value
      if (pop_size <= 0.0) throw std::invalid_argument("population size <= 0");
      pop_size = 1.0 / (2 * pop_size);
    }
    pop_sizes_list_[position].push_back(pop_size);
  }
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
void Model::addPopulationSizes(const double time, const double pop_size,
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
void Model::addPopulationSize(const double time, const size_t pop, double population_size,
                              const bool &time_scaled, const bool &relative) {
  checkPopulation(pop);
  size_t position = addChangeTime(time, time_scaled);
  if (relative) population_size *= default_pop_size();

  if (population_size <= 0.0) throw std::invalid_argument("population size <= 0");
  if (pop_sizes_list_.at(position).empty()) addPopulationSizes(time, nan("value to replace"), time_scaled);
  pop_sizes_list_.at(position).at(pop) = 1.0/(2*population_size);
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
void Model::addGrowthRates(const double time, const std::vector<double> &growth_rates,
                           const bool &time_scaled, const bool &rate_scaled) {
  if ( growth_rates.size() != population_number() )
    throw std::logic_error("Growth rates values do not meet the number of populations");
  size_t position = addChangeTime(time, time_scaled);

  growth_rates_list_[position].clear();
  for (double rate : growth_rates) {
    if (rate_scaled) rate *= scaling_factor();
    growth_rates_list_[position].push_back(rate);
  }
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
void Model::addGrowthRates(const double time, const double growth_rate,
                           const bool &time_scaled, const bool &rate_scaled) {
  addGrowthRates(time, std::vector<double>(population_number(), growth_rate), time_scaled, rate_scaled);
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
void Model::addGrowthRate(const double time, const size_t population,
                          double growth_rate, const bool &time_scaled, const bool &rate_scaled) {
  checkPopulation(population);
  size_t position = addChangeTime(time, time_scaled);
  if (rate_scaled) growth_rate *= scaling_factor();
  if (growth_rates_list_.at(position).empty()) addGrowthRates(time, nan("number to replace"), time_scaled);
  growth_rates_list_.at(position).at(population) = growth_rate;
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
  checkPopulation(source);
  checkPopulation(sink);
  size_t position = addChangeTime(time, scaled_time);
  if (scaled_rates) mig_rate *= scaling_factor();
  if (mig_rates_list_.at(position).empty()) {
    addSymmetricMigration(time, nan("value to replace"), scaled_time);
  }
  mig_rates_list_.at(position).at(getMigMatrixIndex(source, sink)) = mig_rate;
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
  if (scaled_rates) scaling = scaling_factor();
  if ( mig_rates.size() != population_number()*population_number() )
    throw std::logic_error("Migration rates values do not meet the number of populations");

  size_t position = addChangeTime(time, scaled_time);
  mig_rates_list_[position].clear();
  mig_rates_list_[position].reserve(popnr*popnr-popnr);
  for (size_t i = 0; i < popnr; ++i) {
    for (size_t j = 0; j < popnr; ++j) {
      if (i == j) continue;
      mig_rates_list_[position].push_back(mig_rates.at(i*popnr+j) * scaling);
    }
  }
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
void Model::addSymmetricMigration(const double time, const double mig_rate,
                                  const bool &time_scaled, const bool &rate_scaled) {
    std::vector<double> mig_rates = std::vector<double>(population_number()*population_number(), mig_rate);
    this->addMigrationRates(time, mig_rates, time_scaled, rate_scaled);
  }


void Model::addSingleMigrationEvent(const double time, const size_t source_pop,
                                    const size_t sink_pop, const double fraction,
                                    const bool &time_scaled) {

  size_t position = addChangeTime(time, time_scaled);
  size_t popnr = population_number();

  if ( time < 0.0 ) throw std::invalid_argument("Single migration event: Negative time");
  if ( source_pop >= population_number() ) throw std::invalid_argument("Single migration event: Unknown population");
  if ( sink_pop >= population_number() ) throw std::invalid_argument("Single migration event: Unknown population");
  if ( fraction < 0.0 || fraction > 1.0 ) throw std::invalid_argument("Single migration event: Fraction out of range");

  if ( single_mig_probs_list_.at(position).empty() ) {
    single_mig_probs_list_.at(position) = std::vector<double>(popnr*popnr-popnr, 0.0);
  }

  single_mig_probs_list_.at(position).at(getMigMatrixIndex(source_pop, sink_pop)) = fraction;
  this->has_migration_ = true;
}


std::ostream& operator<<(std::ostream& os, Model& model) {
  size_t n_pops = model.population_number();
  os << "---- Model: ------------------------" << std::endl;
  os << "Total Sample Size: " << model.sample_size() << std::endl;
  os << "N0 is assumed to be " << model.default_pop_size() << std::endl;

  model.resetSequencePosition();
  for (size_t idx = 0; idx < model.countChangePositions(); ++idx) {
    os << std::endl << "At Position " << model.getCurrentSequencePosition() << ":" << std::endl;
    os << " Mutation Rate: " << model.mutation_rate() << std::endl;
    os << " Recombination Rate: " << model.recombination_rate() << std::endl;
    model.increaseSequencePosition();
  }
  model.resetSequencePosition();

  model.resetTime();
  for (size_t idx = 0; idx < model.countChangeTimes(); ++idx) {
    os << std::endl << "At Time " << model.getCurrentTime() << ":" << std::endl;

    os << " Population Sizes: ";
    for (size_t pop = 0; pop < n_pops; ++pop)
      os << std::setw(10) << std::right << model.population_size(pop, model.getCurrentTime());
    os << std::endl;

    os << " Growth Rates:     ";
    for (size_t pop = 0; pop < n_pops; ++pop)
      os << std::setw(10) << std::right << model.growth_rate(pop);
    os << std::endl;

    os << " Migration Matrix: " << std::endl;
    for (size_t i = 0; i < n_pops; ++i) {
      for (size_t j = 0; j < n_pops; ++j) {
        os << std::setw(10) << std::right << model.migration_rate(i, j);
      }
      os << std::endl;
    }

    for (size_t i = 0; i < n_pops; ++i) {
      for (size_t j = 0; j < n_pops; ++j) {
        if (model.single_mig_pop(i, j) != 0) {
          os << " " << model.single_mig_pop(i, j) * 100 << "% of pop "
             << i + 1 << " move to pop " << j + 1 << std::endl;
        }
      }
    }

    if (idx < model.countChangeTimes() - 1) model.increaseTime();
  }
  model.resetTime();
  os << "------------------------------------" << std::endl;
  return(os);
}


void Model::updateTotalMigRates(const size_t position) {
  if ( total_mig_rates_list_.at(position).empty() ) {
    total_mig_rates_list_.at(position) = std::vector<double>(population_number(), 0.0);
  }

  std::vector<double>* mig_rates = &(total_mig_rates_list_.at(position));

  for (size_t i = 0; i < population_number(); ++i) {
    for (size_t j = 0; j < population_number(); ++j) {
      if (i == j) continue;
      mig_rates->at(i) += mig_rates_list_.at(position).at( getMigMatrixIndex(i,j) );
    }
    if (mig_rates->at(i) > 0) has_migration_ = true;
  }
}


void Model::finalize() {
  fillVectorList(mig_rates_list_, default_mig_rate);
  fillVectorList(growth_rates_list_, default_growth_rate);
  calcPopSizes();

  for (size_t j = 0; j < mig_rates_list_.size(); ++j) {
    if (mig_rates_list_.at(j).empty()) continue;
    updateTotalMigRates(j);
  }

  // Fill in missing recombination rates
  assert( mutation_rates_.at(0) != -1 );
  assert( recombination_rates_.at(0) != -1 );
  for (size_t j = 1; j < change_position_.size(); ++j) {
    if (mutation_rates_.at(j) == -1) {
      mutation_rates_.at(j) = mutation_rates_.at(j-1);
    }

    if (recombination_rates_.at(j) == -1) {
      recombination_rates_.at(j) = recombination_rates_.at(j-1);
    }
  }

  resetTime();
  resetSequencePosition();
  check();
}


void Model::calcPopSizes() {
  // Set initial population sizes
  if (pop_sizes_list_.at(0).empty()) addPopulationSizes(0, default_pop_size());
  else {
    // Replace values not set by the user with the default size
    for (size_t pop = 0; pop < population_number(); ++pop) {
      if (std::isnan(pop_sizes_list_.at(0).at(pop)))
        addPopulationSize(0, pop, default_pop_size());
    }
  }

  size_t last_growth = -1;
  double duration = -1;
  for (size_t i = 1; i < change_times_.size(); ++i) {
    if (! growth_rates_list_.at(i - 1).empty()) last_growth = i - 1;

    // Make sure we always have a pop sizes vector
    if (pop_sizes_list_.at(i).empty()) {
      addPopulationSizes(change_times_.at(i), nan("value to replace"));
      assert( ! pop_sizes_list_.at(i).empty() );
    }

    // Calculate the effective duration of a growth period if it ends here
    duration = change_times_.at(i) - change_times_.at(i - 1);

    // Calculate the population sizes:
    for (size_t pop = 0; pop < population_number(); ++pop) {
      // If the user explicitly gave a value => always use this value
      if ( !std::isnan(pop_sizes_list_.at(i).at(pop)) ) continue;

      assert(!std::isnan(pop_sizes_list_.at(i - 1).at(pop)));
      // Otherwise use last value
      pop_sizes_list_.at(i).at(pop) = pop_sizes_list_.at(i - 1).at(pop);

      // ... and scale it if there was growth
      if (last_growth != -1) {
        pop_sizes_list_.at(i).at(pop) *=
          std::exp((growth_rates_list_.at(last_growth).at(pop)) * duration);
      }
    }
  }
}


void Model::check() {
  // Sufficient sample size?
  if (sample_size() < 2) throw std::invalid_argument("Sample size needs be to at least 2");

  // Structure without migration?
  if (population_number() > 1 && !has_migration())
    throw std::invalid_argument("Model has multiple populations but no migration. Coalescence impossible");

  // Are the bias heights in the correct order?
  for( size_t idx = 1; idx < bias_heights().size(); idx++ ){
    if ( bias_heights()[idx] < bias_heights()[idx-1]  ) {
      throw std::invalid_argument(std::string("The bias heights must be input in order, recent to ancient"));
    }
  }

  // Are bias heights and bias strengths compatible?
  if ( bias_heights().size() != bias_strengths().size() + 1 ) {
    throw std::invalid_argument(std::string("the input bias_strengths should have one more value than bias_heights") +
                                std::string(" as bias_heights declares the time boundaries between bias_strengths"));
  }
}


void Model::fillVectorList(std::vector<std::vector<double> > &vector_list, const double default_value) {
  std::vector<double>* last = NULL;
  std::vector<double>* current = NULL;
  for (size_t j = 0; j < vector_list.size(); ++j) {
    current = &(vector_list.at(j));
    if (current->empty()) continue;

    for (size_t i = 0; i < current->size(); ++i) {
      if ( !std::isnan(current->at(i)) ) continue;

      if (last == NULL) current->at(i) = default_value;
      else current->at(i) = last->at(i);
    }
    last = current;
  }
}


void Model::addPopulation() {
  // Create the new population
  size_t new_pop = population_number();
  this->set_population_number(new_pop+1);

  // Change Vectors
  addPopToVectorList(growth_rates_list_);
  addPopToVectorList(pop_sizes_list_);

  // Change Matrices
  addPopToMatrixList(mig_rates_list_, new_pop);
  addPopToMatrixList(single_mig_probs_list_, new_pop, 0);
}


void Model::addPopToMatrixList(std::vector<std::vector<double> > &vector_list, size_t new_pop, double default_value) {
  for (auto it = vector_list.begin(); it!= vector_list.end(); ++it) {
    if (it->empty()) continue;
    for (size_t i = 0; i < new_pop; ++i) {
      it->insert(it->begin() + getMigMatrixIndex(i, new_pop), default_value);
    }
    for (size_t i = 0; i < new_pop; ++i) {
      it->insert(it->begin() + getMigMatrixIndex(new_pop, i), default_value);
    }
  }
}


void Model::addPopToVectorList(std::vector<std::vector<double> > &vector_list) {
  for (auto it = vector_list.begin(); it!= vector_list.end(); ++it) {
    if (it->empty()) continue;
    it->push_back(nan("value to replace"));
  }
}


/**
 * @brief Sets the recombination rate
 *
 * @param rate The recombination rate per sequence length per time unit.
 * The sequence length can be either per locus or per base pair and the time
 * can be given in units of 4N0 generations ("scaled") or per generation.
 *
 * @param loci_length The length of every loci.
 * @param per_locus Set to TRUE, if the rate is given per_locus, and to FALSE
 * if the rate is per base pair.
 * @param scaled Set to TRUE is the rate is scaled with 4N0, or to FALSE if
 * it isn't
 */
void Model::setRecombinationRate(double rate,
                                 const bool &per_locus,
                                 const bool &scaled,
                                 const double seq_position) {

  if (rate < 0.0) {
    throw std::invalid_argument("Recombination rate must be non-negative");
  }

  if (scaled) rate /= 4.0 * default_pop_size();
  if (per_locus) {
    if (loci_length() <= 1) {
      throw std::invalid_argument("Locus length must be at least two");
    }
    rate /= loci_length()-1;
  }

  if (rate > 0.0) has_recombination_ = true;
  recombination_rates_[addChangePosition(seq_position)] = rate;
}


/**
 * @brief Sets the mutation rate
 *
 * @param rate The mutation rate per locus, either scaled as theta=4N0*u or
 * unscaled as u.
 * @param per_locus TRUE if the rate is per locus, FALSE if per base.
 * @param scaled Set this to TRUE if you give the mutation rate in scaled
 * units (e.g. as theta rather than as u).
 */
void Model::setMutationRate(double rate, const bool &per_locus, const bool &scaled,
                            const double seq_position) {
  if (scaled) rate /= 4.0 * default_pop_size();

  size_t idx = addChangePosition(seq_position);
  if (per_locus) {
    mutation_rates_.at(idx) = rate / loci_length();
  } else {
    mutation_rates_.at(idx) = rate;
  }
}


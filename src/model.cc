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

//Model::Model() : default_growth_rate(0.0), default_mig_rate(0.0)
//{ 
  //this->init();
//}

Model::Model(size_t sample_size) {
  this->init();

  this->addSampleSizes(0.0, std::vector<size_t>(1, sample_size));
  this->resetTime();
}

//Model::Model(size_t sample_size) : default_growth_rate(0.0), default_mig_rate(0.0)
//{
  //this->init();

  //this->addSampleSizes(0.0, std::vector<size_t>(1, sample_size));
  //this->resetTime();
//}

Model::~Model() { 
  //std::cout << "Called ~Model" << std::endl;
  deleteParList(pop_sizes_list_);
  deleteParList(growth_rates_list_);
  deleteParList(mig_rates_list_);
  deleteParList(total_mig_rates_list_);
  deleteParList(single_mig_probs_list_);
}

void Model::init() {
  this->addChangeTime(0.0);

  this->set_population_number(1);
  this->set_loci_number(1);
  this->set_loci_length(100000);
  this->set_mutation_rate(0.0);
  this->set_recombination_rate(0.0);

  this->set_exact_window_length(0);
  this->set_prune_interval(0);
  this->set_mutation_exact_number(-1); //-1 is equivalent to infinity

  this->resetTime();
}


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


void Model::addPopulationSizes(double time, const std::vector<double> &pop_sizes, bool relative) {
  if ( pop_sizes.size() != population_number() ) 
    throw std::logic_error("Population size values do not meet the number of populations");
  auto pop_sizes_heap = new std::vector<double>(pop_sizes);
  if (relative) {
    for (auto it = pop_sizes_heap->begin(); it != pop_sizes_heap->end(); ++it) {
      if (isnan(*it)) continue;
      else *it *= this->default_pop_size; 
    }
  }
  size_t position = addChangeTime(time);
  pop_sizes_list_[position] = pop_sizes_heap;  
}


void Model::addPopulationSizes(const double &time, const double &pop_size, bool relative) {
  addPopulationSizes(time, std::vector<double>(population_number(), pop_size), relative);
}
  

void Model::addPopulationSize(const double &time, const size_t &pop, const double &population_size, bool relative) {
  size_t position = addChangeTime(time);
  if (pop_sizes_list_.at(position) == NULL) addPopulationSizes(time, nan("value to replace"));
  pop_sizes_list_.at(position)->at(pop) = population_size*default_pop_size;
}


void Model::addGrowthRates(const double &time, const std::vector<double> &growth_rates) {
  if ( growth_rates.size() != population_number() ) 
    throw std::logic_error("Growth rates values do not meet the number of populations");
  std::vector<double>* growth_rates_heap = new std::vector<double>(growth_rates);
  size_t position = addChangeTime(time);
  growth_rates_list_[position] = growth_rates_heap; 

}

void Model::addGrowthRates(const double &time, const double &growth_rate) {
  addGrowthRates(time, std::vector<double>(population_number(), growth_rate));
}

void Model::addGrowthRate(const double &time, const size_t &population, const double &growth_rate) {
  size_t position = addChangeTime(time);
  if (growth_rates_list_.at(position) == NULL) addGrowthRates(time, nan("number to replace")); 
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
 */
void Model::addMigrationRate(double time, size_t source, size_t sink, double mig_rate) {
  size_t position = addChangeTime(time);
  if (mig_rates_list_.at(position) == NULL) addSymmetricMigration(time, nan("value to replace")); 
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
 */
void Model::addMigrationRates(double time, const std::vector<double> &mig_rates) {
  double popnr = population_number();
  if ( mig_rates.size() != population_number()*population_number() ) 
    throw std::logic_error("Migration rates values do not meet the number of populations");
  std::vector<double>* mig_rates_heap = new std::vector<double>();
  mig_rates_heap->reserve(popnr*popnr-popnr);
  for (size_t i = 0; i < popnr; ++i) {
    for (size_t j = 0; j < popnr; ++j) {
      if (i == j) continue;
      mig_rates_heap->push_back(mig_rates.at(i*popnr+j)); 
    }
  }

  size_t position = addChangeTime(time);
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
 */
void Model::addSymmetricMigration(const double &time, const double &mig_rate) {
  std::vector<double> mig_rates = std::vector<double>(population_number()*population_number(), mig_rate);
  this->addMigrationRates(time, mig_rates);
}


void Model::addSingleMigrationEvent(const double &time, const size_t &source_pop, 
                                    const size_t &sink_pop, const double &fraction) {
  
  size_t position = addChangeTime(time);
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
}

void Model::fillVectorList(std::vector<std::vector<double>*> &vector_list, const double &default_value) {
  std::vector<double>* last = NULL; 
  std::vector<double>* current = NULL; 
  for (size_t j = 0; j < vector_list.size(); ++j) {
    current = vector_list.at(j);
    if (current == NULL) continue;

    for (size_t i = 0; i < current->size(); ++i) {
      if ( !isnan(current->at(i)) ) continue;

      if (last == NULL) (current)->at(i) = default_value;
      else current->at(i) = last->at(i); 
    }
    last = current;
  }
}

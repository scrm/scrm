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


void Model::addPopulationSizes(double time, const std::vector<size_t> &pop_sizes) {
  if ( pop_sizes.size() != population_number() ) 
    throw std::logic_error("Population size values do not meet the number of populations");
  std::vector<size_t>* pop_sizes_heap = new std::vector<size_t>(pop_sizes);
  size_t position = addChangeTime(time);
  pop_sizes_list_[position] = pop_sizes_heap;  
}


void Model::addRelativePopulationSizes(double time, const std::vector<double> &population_sizes) {
  std::vector<size_t> abs_pop_sizes;
  std::vector<double>::const_iterator it;
  for (it = population_sizes.begin(); it != population_sizes.end(); ++it) {
    abs_pop_sizes.push_back( *it * this->default_pop_size ); 
  }
  this->addPopulationSizes(time, abs_pop_sizes);
}


void Model::addGrowthRates(double time, const std::vector<double> &growth_rates) {
  if ( growth_rates.size() != population_number() ) 
    throw std::logic_error("Growth rates values do not meet the number of populations");
  std::vector<double>* growth_rates_heap = new std::vector<double>(growth_rates);
  size_t position = addChangeTime(time);
  growth_rates_list_[position] = growth_rates_heap; 

}

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
  updateTotalMigRates(position);
}


std::ostream& operator<<(std::ostream& os, const Model& model) {
  os << "---- Model: ------------------------" << std::endl;
  os << "Mutation rate: " << model.mutation_rate() << std::endl;  
  os << "Recombination rate: " << model.recombination_rate() << std::endl;  
  os << model.pop_sizes_list_ << std::endl << model.growth_rates_list_ << std::endl;
  
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
      mig_rates->at(j) += mig_rates_list_.at(position)->at( getMigMatrixIndex(i,j) );
    }
  }
}; 

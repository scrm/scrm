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
}

void Model::init() {
  this->addChangeTime(0.0);

  this->set_population_number(1);
  this->set_loci_number(1);
  this->set_mutation_rate(0.0);
  this->set_recombination_rate(0.0);

  this->set_exact_window_length(0);
  this->set_prune_interval(0);

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


void Model::addGrowthRates(double time, const std::vector<double> &growth_rates) {
  if ( growth_rates.size() != population_number() ) 
    throw std::logic_error("Growth rates values do not meet the number of populations");
  std::vector<double>* growth_rates_heap = new std::vector<double>(growth_rates);
  size_t position = addChangeTime(time);
  growth_rates_list_[position] = growth_rates_heap; 
}


void Model::print(std::ostream &os) const {
  os << "---- Model: ------------------------" << std::endl;
  os << "Mutation rate: " << this->mutation_rate() << std::endl;  
  os << "Recombination rate: " << this->recombination_rate() << std::endl;  
  for (size_t idx = 0; idx < change_times_.size(); ++idx) { 
    os << std::endl << "At time " << change_times_.at(idx) << ":" << std::endl;  
    if (pop_sizes_list_.at(idx) != NULL) {
      os << " Population sizes: "; 
      printVector(*(pop_sizes_list_.at(idx)), os);
      os << std::endl;
    }
    if (growth_rates_list_.at(idx) != NULL) {
      os << " Growth Rate: "; 
      printVector(*(growth_rates_list_.at(idx)), os);
      os << std::endl;
    }
  }
  os << "------------------------------------" << std::endl;
}

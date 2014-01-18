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


#include "param.h"

void Param::init(){
  this->set_random_seed(-1);
  //this->rho=0;	//double, in ms, rho = 4 * npop * recomb_rate_persite * (nsites-1)
  //this->recomb_rate_persite=0;
  this->log_bool = false;
  this->tree_bool = false;
  this->log_NAME = "scrm.log";
  this->tree_NAME = "scrm.treefile";
  this->set_seg_bool(true);
  this->tmrca_bool = false;
  this->tmrca_NAME = "scrm.tmrcafile";
}


std::ostream& operator<< (std::ostream& stream, const Param& param) {
  //stream << param.argv_[0];
  for (int i = 0; i < param.argc_; ++i) {
    stream << " " << param.argv_[i];
  }
  return stream;
}


void Param::nextArg(std::string option) {
  ++argc_i;
  if (argc_i >= argc_) {
    throw std::invalid_argument(std::string("Not enough parameters when parsing option ") + option);
  }
}

/*! 
 * \brief Read in ms parameters and convert to scrm parameters
 * The first parameters followed by -eG, -eg, -eN, -en, -em, -ema, -es 
 * and -ej options in ms are time t in unit of 4N_0 generations. 
 * In scrm, we define time t in number of generations. 
 */

void Param::parse(Model &model) {
  model = Model();
  size_t sample_size = 0;
  double time = 0;

  if ( argc_ == 0 ) return;
  if (directly_called_){
    argc_i=0;
  }
  else{
    argc_i=1;
    std::cout<<"Indirectely called"<<std::endl;
  }

  this->init();

  while( argc_i < argc_ ){
    //std::cout << argv_[argc_i] << std::endl;
    std::string argv_i = argv_[argc_i];

    if (argc_i == 0) {
      nextArg("nsam (1st option)");
      sample_size = readInput<size_t>(argv_[argc_i]);
      nextArg("howmany (2nd option)");
      model.set_loci_number(readInput<int>(argv_[argc_i]));
    }

    // ------------------------------------------------------------------
    // Mutation 
    // ------------------------------------------------------------------
    else if (argv_i == "-t") {
      nextArg(argv_i);
      model.set_mutation_rate(readInput<double>(argv_[argc_i]), true, true);
    }

    else if (argv_i == "-s"){
      nextArg(argv_i);
      model.set_mutation_exact_number(readInput<size_t>(argv_[argc_i]));
    }

    // ------------------------------------------------------------------
    // Recombination 
    // ------------------------------------------------------------------
    else if (argv_i == "-r") {
      double rec_rate = readNextInput<double>();
      size_t loci_length = readNextInput<size_t>();
      model.set_recombination_rate(rec_rate, loci_length, true, true);
    }

    // ------------------------------------------------------------------
    // Subpopulations 
    // ------------------------------------------------------------------
    // Set number of subpopulations and samples at time 0
    else if (argv_i == "-I") {
      nextArg(argv_i);
      model.set_population_number(readInput<size_t>(argv_[argc_i]));
      std::vector<size_t> sample_size;
      for (size_t i = 0; i < model.population_number(); ++i) {
        nextArg(argv_i);
        sample_size.push_back(readInput<size_t>(argv_[argc_i]));
      }
      model.addSampleSizes(0.0, sample_size);
    }

    // Add samples at arbitrary times
    else if (argv_i == "-eI") {
      nextArg(argv_i);
      time = readInput<double>(argv_[argc_i]);
      std::vector<size_t> sample_size;
      for (size_t i = 0; i < model.population_number(); ++i) {
        nextArg(argv_i);
        sample_size.push_back(readInput<size_t>(argv_[argc_i]));
      }
      model.addSampleSizes(time, sample_size);
    }

    // ------------------------------------------------------------------
    // Populations sizes 
    // ------------------------------------------------------------------
    else if (argv_i == "-eN" || argv_i == "-N") {
      if (argv_i == "-eN") time = readNextInput<double>();
      else time = 0.0;
      model.addPopulationSizes(time, readNextInput<double>(), true, true); 
      model.addGrowthRates(time, 0.0, true);
    }

    else if (argv_i == "-en" || argv_i == "-n") {
      if (argv_i == "-en") time = readNextInput<double>();
      else time = 0.0;
      size_t pop = readNextInput<size_t>() - 1;
      model.addPopulationSize(time, pop, readNextInput<double>(), true, true);
      model.addGrowthRate(time, pop, 0.0, true);
    }

    // ------------------------------------------------------------------
    // Exponential Growth 
    // ------------------------------------------------------------------
    else if (argv_i == "-G" || argv_i == "-eG") {
      if (argv_i == "-eG") time = readNextInput<double>();
      else time = 0.0;
      model.addGrowthRates(time, readNextInput<double>(), true); 
    }

    else if (argv_i == "-g" || argv_i == "-eg") {
      if (argv_i == "-eg") time = readNextInput<double>();
      else time = 0.0;
      size_t pop = readNextInput<size_t>() - 1;
      model.addGrowthRate(time, pop, readNextInput<double>(), true); 
    }

    // ------------------------------------------------------------------
    // Migration
    // ------------------------------------------------------------------
    else if (argv_i == "-ma" || argv_i == "-ema") {
      if (argv_i == "-ema") {
        nextArg(argv_i);
        time = readInput<double>(argv_[argc_i]);
      }
      else time = 0.0;
      std::vector<double> migration_rates;
      for (size_t i = 0; i < model.population_number(); ++i) {
        for (size_t j = 0; j < model.population_number(); ++j) {
          nextArg(argv_i);
          if (i==j) migration_rates.push_back(0.0);
          else migration_rates.push_back(readInput<double>(argv_[argc_i]));
        }
      }
      model.addMigrationRates(time, migration_rates, true, true);
    }

    else if (argv_i == "-M" || argv_i == "-eM") {
      if (argv_i == "-eM") {
        nextArg(argv_i);
        time = readInput<double>(argv_[argc_i]);
      }
      else time = 0.0;
      nextArg(argv_i);
      model.addSymmetricMigration(time, readInput<double>(argv_[argc_i])/(model.population_number()-1), true, true);
    }

    else if (argv_i == "-esme") {
      time = readNextInput<double>();
      size_t source_pop = readNextInput<size_t>() - 1;
      size_t sink_pop = readNextInput<size_t>() - 1;
      double fraction = readNextInput<double>();
      model.addSingleMigrationEvent(time, source_pop, sink_pop, fraction, true); 
    }

    else if (argv_i == "-ej") {
      time = readNextInput<double>();
      size_t source_pop = readNextInput<size_t>() - 1;
      size_t sink_pop = readNextInput<size_t>() - 1;
      model.addSingleMigrationEvent(time, source_pop, sink_pop, 1.0, true); 
      for (size_t i = 0; i < model.population_number(); ++i) {
        if (i == source_pop) continue;
        model.addMigrationRate(time, i, source_pop, 0.0, true);
      }
    }

    // ------------------------------------------------------------------
    // Pruning 
    // ------------------------------------------------------------------
    else if (argv_i == "-l"){
      nextArg(argv_i);
      model.set_exact_window_length(readInput<size_t>(argv_[argc_i]));
    }

    // ------------------------------------------------------------------
    // Output 
    // ------------------------------------------------------------------
    else if (argv_i == "-T"){
      //treefile = argv_[++argc_i];
      tree_bool = true;
    }

    else if (argv_i == "-L"){
      //treefile = argv_[++argc_i];
      tmrca_bool = true;
    }

    else if (argv_i == "-seed"){
      nextArg(argv_i);
      this->set_random_seed( readInput<size_t>(argv_[argc_i]) );
    }

    else if (directly_called_){
      throw std::invalid_argument(std::string("unknown/unexpected argument: ") + argv_i);
    }

    ++argc_i;
  }

  if (model.sample_size() == 0) {
    model.addSampleSizes(0.0, std::vector<size_t>(1, sample_size));
  } 
  else if (model.sample_size() != sample_size) {
    throw std::invalid_argument("Sum of samples not equal to the total sample size");
  }

  if (directly_called_){
    model.finalize();
  }
  model.resetTime();
}



void Param::print_param(){
}	

void Param::log_param(){
  //std::ofstream log_file;
  //log_file.open (log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
  std::ofstream log_file(this->log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
  log_file<<"scrm parameters: \n";
  log_file<<std::setw(10)<<"seed:"<<" "<<std::setw(10)<<random_seed()<< "\n";
  log_file.close();
}		


void print_options(){
  std::cout << "  Options (ms-compatible):" << std::endl;
  std::cout<<std::setw(20)<<"-t theta          Set theta = 4*N0*mu, mu = locus mutation rate per generation" << std::endl;
  std::cout<<std::setw(20)<<"-r rho nsites     Set rho = 4*N0*r, r = locus recombination rate; and locus length" << std::endl;
  std::cout<<std::setw(20)<<"-T                Output gene trees" << std::endl;
  std::cout<<std::setw(20)<<"-s segsites       Sample a fixed number of segregating sites onto the sampled ARG" << std::endl;
  std::cout<<std::setw(20)<<"-I n n1 .. nn     Define n populations, each of effective size N0, each having n1,..,nn samples. (Migration not yet implemented.)" << std::endl;
  std::cout<<std::setw(20)<<"-eN t size        Modify population sizes, to become size*N0 prior to time t" << std::endl;
  std::cout<<std::endl;
  std::cout<<"  Further options:" << std::endl;
  std::cout<<std::setw(20)<<"-seed SEED        Set seed for random number generator (default: use system time)"<<std::endl;
  std::cout<<std::setw(20)<<"-eI t n n1 .. nn  As -I but define  n  new populations at time t in past" << std::endl;
  //std::cout<<std::setw(20)<<"-l exact_window_length"<<"  --  "
  //    <<"User define the length of the exact window."<<std::endl;
  //std::cout<<std::setw(20)<<"-log [LOGFILE]"<<"  --  "<< "User specify the log file name, scrm.log by default."<<std::endl;
  //std::cout<<std::setw(20)<<"-T myTREEFILE"<<"  --  "<< "User specify the tree file name, TREEFILE by default."<<std::endl;
  std::cout<<std::endl;
}


void print_help(){
  std::cout << "usage: scrm nsam howmany" << std::endl;
  print_options();
  print_example();
}


void print_example(){	
  std::cout<<"Examples:"<<std::endl;
  std::cout<<"./scrm 3 1"<<std::endl;
  std::cout<<"./scrm 6 3 -t 10 -r 40 -npop 20000 "<<std::endl;
  std::cout<<"./scrm 5 3 -t 20 -r 30 -npop 10000 -seed 1314 -log"<<std::endl;
  std::cout<<"./scrm 6 1 -t 30 -r 40 10000 -nsites 2000 -log"<<std::endl;
  std::cout<<"./scrm 6 2 -t 25 -r 30 -nsites 10000 -log LOGFILE"<<std::endl;
  std::cout<<"./scrm 6 2 -t 40 -r 40 -log LOGFILE -T mytree"<<std::endl;
}

void appending_log_file(std::string log_file_NAME,std::string log_file_input /*! Information added*/){
  std::ofstream log_file;
  log_file.open (log_file_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
  log_file << log_file_input << "\n";
  log_file.close();
}

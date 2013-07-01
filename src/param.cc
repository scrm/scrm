/* 
 * param.cc
 */

#include "param.h"

void Param::init(){
  this->random_seed=time(0);
  //this->rho=0;	//double, in ms, rho = 4 * npop * recomb_rate_persite * (nsites-1)
  //this->recomb_rate_persite=0;
  this->log_bool=false;
  this->log_NAME="scrm.log";
  this->treefile="TREEFILE";
  this->set_seg_bool(true);
  this->total_mut=0;
  this->set_tmrca_bool(false);
  this->tmrca_NAME="tmrcaFILE";
}

Model* Param::parse() {
  Model* model = new Model();
  size_t sample_size = 0;
  double time = 0;
  if ( argc_ == 0 ) return model;

  int argc_i=0;
  this->init();

  while( argc_i < argc_ ){
    //std::cout << argv_[argc_i] << std::endl;
    std::string argv_i = argv_[argc_i];

    if (argc_i == 0) {
      sample_size = readInput<size_t>(argv_[++argc_i]);
      model->set_loci_number(readInput<int>(argv_[++argc_i]));
    }	

    else if (argv_i == "-t") {
      model->set_mutation_rate(readInput<double>(argv_[++argc_i]));
    }

    else if (argv_i == "-r") {
      model->set_recombination_rate(readInput<double>(argv_[++argc_i]));
      model->set_loci_length(readInput<size_t>(argv_[++argc_i]));
    }

    else if (argv_i == "-I") {
      model->set_population_number(readInput<size_t>(argv_[++argc_i]));
      std::vector<size_t> sample_size;
      for (size_t i = 0; i < model->population_number(); ++i) {
        sample_size.push_back(readInput<size_t>(argv_[++argc_i]));
      }
      model->addSampleSizes(0.0, sample_size);
    }

    else if (argv_i == "-eI") {
      time = readInput<double>(argv_[++argc_i]);
      std::vector<size_t> sample_size;
      for (size_t i = 0; i < model->population_number(); ++i) {
        sample_size.push_back(readInput<size_t>(argv_[++argc_i]));
      }
      model->addSampleSizes(time, sample_size);
    }

    else if (argv_i == "-eN") {
      time = readInput<double>(argv_[++argc_i]);
      std::vector<double> pop_sizes(model->population_number(),
                                    readInput<double>(argv_[++argc_i]));
      model->addRelativePopulationSizes(time, pop_sizes); 
      //std::cout << model << std::endl;
    }

    else if (argv_i == "-l"){
      model->set_exact_window_length(readInput<size_t>(argv_[++argc_i]));
    }

    else if (argv_i == "-seed"){
      random_seed = readInput<int>(argv_[++argc_i]);
    }

    else if (argv_i == "-T"){
      treefile = argv_[++argc_i];
    }

    else {
      throw std::invalid_argument(std::string("unknown/unexpected argument: ") + argv_i);
    }

    ++argc_i; 
    /* 


    if (argv_i=="-s"){
      read_input_to_param<int>(argv[argc_i+1],total_mut);
      argc_i++;
    }

    if (argv_i=="-tmrca"){
      tmrca_bool=true;
      argc_i++;
      if (argc_i < argc){
        if (argv[argc_i][0]!='-'){
          tmrca_NAME=argv[argc_i];
          //argc_i++;
        }
        else{argc_i--;}
      }
    }

    if (argv_i=="-log"){
      log_bool=true;
      argc_i++;
      if (argc_i < argc){
        if (argv[argc_i][0]!='-'){
          log_NAME=argv[argc_i];
          //argc_i++;
        }
        else{argc_i--;}
      }			
    }
    argc_i++;
  */
  }

  if (model->sample_size() == 0) {
    //std::cout << "adding samples" << std::endl;
    model->addSampleSizes(0.0, std::vector<size_t>(1, sample_size));
  } 
  else if (model->sample_size() != sample_size) {
    throw std::invalid_argument("Sum of samples not equal to the total sample size");
  }

  //std::cout << "resetting time" << std::endl;
  model->resetTime();

  /*
  recomb_rate_persite=rho/4 / npop / (nsites-1);
  remove(treefile.c_str());
  if (log_bool){
    remove(log_NAME.c_str());
    log_param();
  }
  if (tmrca_bool){remove(tmrca_NAME.c_str());}
  */
  return model;
}



void Param::print_param(){
}	

void Param::log_param(){
  //std::ofstream log_file;
  //log_file.open (log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
  std::ofstream log_file(this->log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
  log_file<<"scrm parameters: \n";
  log_file<<std::setw(10)<<"seed:"<<" "<<std::setw(10)<<random_seed<< "\n";
  log_file.close();
}		

void print_help(){
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<"*****************************************************************"<<std::endl;
  std::cout<<"*			scrm					*"<<std::endl;
  std::cout<<"*		  Author: Paul R Staab, Sha Zhu			*"<<std::endl;
  std::cout<<"*****************************************************************"<<std::endl;
  std::cout<<std::endl<<std::endl;
  std::cout<<"Too few command line arguments"<<std::endl;
  std::cout<<"usage: ms nsam howmany"<<std::endl;
  std::cout<<"	Options:"<<std::endl;
  print_option();
  print_example();
  exit(1);
}

void print_option(){
  //std::cout<<std::setw(20)<<"-h or -help"<<"  --  "<<"Help. List the following content."<<std::endl;
  std::cout<<std::setw(20)<<"-r RHO"<<"  --  "<<"User define the recombination rate RHO, per gerneration per site."<<std::endl;
  std::cout<<std::setw(20)<<"-nsites NSITES"<<"  --  "<<"User define the sequence length NSITES."<<std::endl;
  std::cout<<std::setw(20)<<"-t THETA"<<"  --  "<<"User define the mutation rate THETA."<<std::endl;
  std::cout<<std::setw(20)<<"-npop NPOP"<<"  --  "<<"User define the population size NPOP."<<std::endl;
  std::cout<<std::setw(20)<<"-seed SEED"<<"  --  "<<"User define the random SEED."<<std::endl;
  std::cout<<std::setw(20)<<"-l exact_window_length"<<"  --  "
      <<"User define the length of the exact window."<<std::endl;
  std::cout<<std::setw(20)<<"-log [LOGFILE]"<<"  --  "<< "User specify the log file name, scrm.log by default."<<std::endl;
  std::cout<<std::setw(20)<<"-T myTREEFILE"<<"  --  "<< "User specify the tree file name, TREEFILE by default."<<std::endl;
}


void print_example(){	
	std::cout<<"Example:"<<std::endl;
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

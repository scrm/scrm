#ifndef SCRM_PARAM_INCLUDED
#define SCRM_PARAM_INCLUDED

#include <iostream>
#include <iomanip>      
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>  
#include <stdio.h>
#include <stdexcept>
#include "mtrand.h"
#include "model.h"

namespace scrm{
  class model;

	class param{
		public:
		bool seg_bool;
		int total_mut;
		int random_seed;
		int nreps; // number of replicates
		double nsites; //number of sites
		
		int npop; // population size Ne, diploid samples
		int nsam; // number of samples
		/*! \todo Check for the definition of rho */ 
		double rho; //Recombination rate: rho = 4*N*r where r is the recombination rate between the ends of the segment being simulated (not the per site rate, per site rate is normally 0.00000001) 
		//and nsites is the number of sites between which recombination can occur.
		double recomb_rate_persite;
		double theta; // mutation rate: theta=4*N*mu, where mutation rate between the two end of the sequence 
		
		size_t ith_change; // the i_th change of the genealogy
		std::string treefile;
		size_t exact_window_length;
		
		bool tmrca_bool;
		std::string tmrca_NAME;
		
		void init();
		param();
		param(Model *model, int argc, char *argv[]);
		
		void print_param();
		
		bool log_bool;
		std::string log_NAME;
		
		void log_param();
		
		
		//private:
		
		
	};
}

namespace scrm{
	//class scrm_help{
	//public:
		void print_help();
		void print_example();
		void print_option();	
	//};
}


//void appending_log_file(std::string log_file_input /*! Information added*/);

template<class T>
void read_input_to_param(char inchar[], T &input)
{
	if (isdigit(inchar[0])){
		std::istringstream para_istrm(inchar);
		para_istrm >> input;
	}
	else{
            throw std::invalid_argument("Invalid argument type. ");
	}	
}

//void read_input_to_double(char inchar[], double &input);
//void read_input_to_int(char inchar[], int &input);
//void read_input_to_size_t(char inchar[], size_t &input);

#endif //SCRM_PARAM_INCLUDED

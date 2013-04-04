/* param.h 
 * */

#include<iostream>
#include<iomanip>      
#include<string>
#include<sstream>
#include<fstream>
#include<iostream>
#include<stdlib.h>  
#include<stdio.h>
class param{
	public:
	
	int random_seed;
	
	int nsites; //number of sites
	
	int npop; // population size Ne
	int nsam; // number of samples
	/*! \todo Check for the definition of rho */ 
	double rho; // recombination rate per site r, recombination rate per generation rho is 4*N*r, where N is the number of population size per generation 
	//Recombination rate: rho = 4*N*r where r is the recombination rate between the ends of the segment being simulated (not the per site rate) 
	//and nsites is the number of sites between which recombination can occur.
	double theta; // mutation rate per site theta, 
	
	size_t ith_change; // the i_th change of the genealogy
	
	param();
	param(int argc, char *argv[]);
	
	void print_param();
	
	bool log_bool;
	std::string log_NAME;
	
	private:
	void log_param();
};


void print_help();
void read_input_to_double(char inchar[], double &input);
void read_input_to_int(char inchar[], int &input);
void read_input_to_size_t(char inchar[], size_t &input);
void appending_log_file(std::string log_file_input /*! Information added*/);


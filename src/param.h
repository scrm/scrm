/* param.h 
 * */

#include<iostream>
#include<iomanip>      


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
	
	param();
	param(int argc, char *argv[]);
	
	void print_param();
};

/* param.cc
 * 
 */
 
#include"param.h"
 
param::param(){
	nsites=25;
	nsam=5;
	npop=10000;
	theta=0.00001;
	rho=0.00002;	
}

param::param(int argc, char *argv[]){
	//size_t arg=0;
     //if( argv[arg][0] != '-' ) { fprintf(stderr," argument should be -%s ?\n", argv[arg]); usage();}

	
	}
	
void param::print_param(){
	std::cout<<std::setw(6)<<"nsites"<<std::setw(10)<<nsites<<std::endl;
	std::cout<<std::setw(6)<<"nsam"<<std::setw(10)<<nsam<<std::endl;
	std::cout<<std::setw(6)<<"npop"<<std::setw(10)<<npop<<std::endl;
	std::cout<<std::setw(6)<<"theta"<<std::setw(10)<<theta<<std::endl;
	std::cout<<std::setw(6)<<"rho"<<std::setw(10)<<rho<<std::endl;
	}	

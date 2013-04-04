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
	ith_change=25;
}

param::param(int argc, char *argv[]){
	int argc_i=0;
	bool help=false;

	random_seed=time(0); //int
	nsites=25; //int
	nsam=5; //int
	npop=10000; //int
	theta=0.00001;	//double
	rho=0.00002;	//double
	ith_change=25; //size_t
	
	while( argc_i < argc ){
		std::string argv_i(argv[argc_i]);
		if (argv_i=="-h" || argv_i=="-help"){
			help=true;
		}
		if (argv_i=="-seed"){
			read_input_to_int(argv[argc_i+1],random_seed);
			argc_i++;
		}
		if (argv_i=="-nsam"){
			read_input_to_int(argv[argc_i+1],nsam);
			argc_i++;
		}	
		if (argv_i=="-npop"){
			read_input_to_int(argv[argc_i+1],npop);
			argc_i++;
		}	
		if (argv_i=="-t"){
			read_input_to_double(argv[argc_i+1],theta);
			argc_i++;
		}
		
		if (argv_i=="-r"){
			read_input_to_double(argv[argc_i+1],rho);
			argc_i++;
		}
		argc_i++;
	}
	
	if (help ){
		print_help();
	}
	remove("scrm.log");
	log_param();
}



void read_input_to_double(char inchar[], double &input){
	if (isdigit(inchar[0])){
		std::istringstream para_istrm(inchar);
		para_istrm >> input;
	}
	else{
		std::cout<<"Error"<<std::endl;//Error message
	}	
}

void read_input_to_int(char inchar[], int &input){
	if (isdigit(inchar[0])){
		std::istringstream para_istrm(inchar);
		para_istrm >> input;
	}
	else{
		std::cout<<"Error"<<std::endl;//Error message
	}	
}

void read_input_to_size_t(char inchar[], size_t &input){
	if (isdigit(inchar[0])){
		std::istringstream para_istrm(inchar);
		para_istrm >> input;
	}
	else{
		std::cout<<"Error"<<std::endl; //Error message
	}	
}
	
	
void param::print_param(){
	std::cout<<std::setw(6)<<"nsites"<<std::setw(10)<<nsites<<std::endl;
	std::cout<<std::setw(6)<<"nsam"<<std::setw(10)<<nsam<<std::endl;
	std::cout<<std::setw(6)<<"npop"<<std::setw(10)<<npop<<std::endl;
	std::cout<<std::setw(6)<<"theta"<<std::setw(10)<<theta<<std::endl;
	std::cout<<std::setw(6)<<"rho"<<std::setw(10)<<rho<<std::endl;
	}	

void param::log_param(){
	std::ofstream log_file;
	log_file.open ("scrm.log", std::ios::out | std::ios::app | std::ios::binary); 
	log_file<<"Simulation parameters: \n";
	log_file<<std::setw(10)<<"nsites ="<<std::setw(10)<<nsites<< "\n";
	log_file<<std::setw(10)<<"nsam ="<<std::setw(10)<<nsam<< "\n";
	log_file<<std::setw(10)<<"npop ="<<std::setw(10)<<npop<< "\n";
	log_file<<std::setw(10)<<"theta ="<<std::setw(10)<<theta<< "\n";
	log_file<<std::setw(10)<<"rho ="<<std::setw(10)<<rho<< "\n";
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
	std::cout<<std::setw(20)<<"-h or -help"<<"  --  "<<"Help. List the following content."<<std::endl;
	std::cout<<std::setw(20)<<"-r RHO"<<"  --  "<<"User define the recombination rate RHO."<<std::endl;
	std::cout<<std::setw(20)<<"-t THETA"<<"  --  "<<"User define the mutation rate THETA."<<std::endl;
	std::cout<<std::setw(20)<<"-nsam NSAM"<<"  --  "<<"User define the sample size NSAM."<<std::endl;
	std::cout<<std::setw(20)<<"-npop NPOP"<<"  --  "<<"User define the population size NPOP."<<std::endl;
	std::cout<<std::setw(20)<<"-seed SEED"<<"  --  "<<"User define the random SEED."<<std::endl;
	std::cout<<"Example:"<<std::endl;
	std::cout<<"./scrm -help"<<std::endl;
	std::cout<<"./scrm -t 0.002 -r 0.00004 -npop 20000 -nsam 6"<<std::endl;
	std::cout<<"./scrm -t 0.0002 -r 0.00003 -npop 10000 -nsam 5 -seed 1314"<<std::endl;
	exit(1);
}


void appending_log_file(std::string log_file_input /*! Information added*/){
	std::ofstream log_file;
	log_file.open ("scrm.log", std::ios::out | std::ios::app | std::ios::binary); 
	log_file << log_file_input << "\n";
	log_file.close();
}



/* param.cc
 * 
 */
 
#include"param.h"
using namespace scrm; 
 
param::param(){
	nsites=25;
	nsam=5;
	npop=10000;
	theta=0.00001;
	rho=0.00000002;	
	ith_change=25;
  exact_window_length=0;
	log_bool=false;
	log_NAME="scrm.log";

}

param::param(int argc, char *argv[]){
	

	
	int argc_i=0;
	//read_input_to_int(argv[argc_i+1],nsam);
	//argc_i++;

	random_seed=time(0); //int
	nsites=25; //int
	//nsam=5; //int
	npop=10000; //int
	theta=0.00001;	//double
	rho=0.00000002;	//double per genration per site, in ms, RHO = 4 * npop * rho * (nsites-1), 
	ith_change=25; //size_t
	exact_window_length=0; //size_t
	log_bool=false;
	log_NAME="scrm.log";
	
	while( argc_i < argc ){
		
		std::string argv_i(argv[argc_i]);
		
		if (argv_i=="-nsam"){ // if scrm is not called, use this option read in the number of samples
			//read_input_to_int(argv[argc_i+1],nsam);
			read_input_to_param<int>(argv[argc_i+1],nsam);
			argc_i++;
		}
		
		if (argv_i=="scrm" || argv_i=="./scrm" || argv_i=="./scrm_dbg" || argv_i=="./scrm_prof"){ // if scrm is directly called
			//read_input_to_int(argv[argc_i+1],nsam);
			read_input_to_param<int>(argv[argc_i+1],nsam);
			argc_i++;
		}
		
		if (argv_i=="-seed"){
			//read_input_to_int(argv[argc_i+1],random_seed);
			read_input_to_param<int>(argv[argc_i+1],random_seed);
			argc_i++;
		}
		
		if (argv_i=="-nsites"){
			//read_input_to_int(argv[argc_i+1],random_seed);
			read_input_to_param<int>(argv[argc_i+1],nsites);
			argc_i++;
		}

		if (argv_i=="-npop"){
			//read_input_to_int(argv[argc_i+1],npop);
			read_input_to_param<int>(argv[argc_i+1],npop);
			argc_i++;
		}	

		if (argv_i=="-t"){
			//read_input_to_double(argv[argc_i+1],theta);
			read_input_to_param<double>(argv[argc_i+1],theta);
			argc_i++;
		}
		
		if (argv_i=="-r"){
			//read_input_to_double(argv[argc_i+1],rho);
			read_input_to_param<double>(argv[argc_i+1],rho);
			argc_i++;
		}
		
		if (argv_i=="-l"){
			//read_input_to_size_t(argv[argc_i+1],exact_window_length);
			read_input_to_param<size_t>(argv[argc_i+1],exact_window_length);
			argc_i++;
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
	}
	
	if (log_bool){
		remove(log_NAME.c_str());
		log_param();
	}
	
}




	
	
void param::print_param(){
	std::cout<<std::setw(6)<<"nsites"<<std::setw(10)<<nsites<<std::endl;
	std::cout<<std::setw(6)<<"nsam"<<std::setw(10)<<nsam<<std::endl;
	std::cout<<std::setw(6)<<"npop"<<std::setw(10)<<npop<<std::endl;
	std::cout<<std::setw(6)<<"theta"<<std::setw(10)<<theta<<std::endl;
	std::cout<<std::setw(6)<<"rho"<<std::setw(10)<<rho<<std::endl;
	std::cout<<std::setw(6)<<"l"<<std::setw(10)<<exact_window_length<<std::endl;
	}	

void param::log_param(){
	std::ofstream log_file;
	log_file.open (log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
	log_file<<"scrm parameters: \n";
	log_file<<std::setw(10)<<"nsites ="<<std::setw(10)<<nsites<< "\n";
	log_file<<std::setw(10)<<"nsam ="<<std::setw(10)<<nsam<< "\n";
	log_file<<std::setw(10)<<"npop ="<<std::setw(10)<<npop<< "\n";
	log_file<<std::setw(10)<<"theta ="<<std::setw(10)<<theta<< "\n";
	log_file<<std::setw(10)<<"rho ="<<std::setw(10)<<rho<< "\n";
	log_file<<std::setw(10)<<"seed:"<<" "<<std::setw(10)<<random_seed<< "\n";
	log_file<<std::setw(10)<<"l ="<<std::setw(10)<<exact_window_length<<std::endl;
	log_file.close();
}		


void scrm::print_help(){
//void scrm_help::print_help(){
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

void scrm::print_option(){
	//std::cout<<std::setw(20)<<"-h or -help"<<"  --  "<<"Help. List the following content."<<std::endl;
	std::cout<<std::setw(20)<<"-r RHO"<<"  --  "<<"User define the recombination rate RHO, per gerneration per site."<<std::endl;
	std::cout<<std::setw(20)<<"-nsites NSITES"<<"  --  "<<"User define the sequence length NSITES."<<std::endl;
	std::cout<<std::setw(20)<<"-t THETA"<<"  --  "<<"User define the mutation rate THETA."<<std::endl;
	std::cout<<std::setw(20)<<"-npop NPOP"<<"  --  "<<"User define the population size NPOP."<<std::endl;
	std::cout<<std::setw(20)<<"-seed SEED"<<"  --  "<<"User define the random SEED."<<std::endl;
	std::cout<<std::setw(20)<<"-l exact_window_length"<<"  --  "
      <<"User define the length of the exact window."<<std::endl;	
}

void scrm::print_example(){	
	std::cout<<"Example:"<<std::endl;
	std::cout<<"./scrm 3"<<std::endl;
	std::cout<<"./scrm 6 -t 0.002 -r 0.00004 -npop 20000 "<<std::endl;
	std::cout<<"./scrm 5 -t 0.0002 -r 0.00003 -npop 10000 -seed 1314"<<std::endl;
	std::cout<<"./scrm 6 -t 0.002 -log -r 0.00004 "<<std::endl;
	std::cout<<"./scrm 6 -t 0.002 -r 0.00004 -log LOGFILE"<<std::endl;
}

void appending_log_file(std::string log_file_NAME,std::string log_file_input /*! Information added*/){
	std::ofstream log_file;
	log_file.open (log_file_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
	log_file << log_file_input << "\n";
	log_file.close();
}
// The rest of the code can be removed ...

void read_input_to_double(char inchar[], double &input){
	if (isdigit(inchar[0])){
		std::istringstream para_istrm(inchar);
		para_istrm >> input;
	}
	else{
		throw std::invalid_argument("Invalid argument type. double");
//		std::cout<<"Error"<<std::endl;//Error message
	}	
}

void read_input_to_int(char inchar[], int &input){
	if (isdigit(inchar[0])){
		std::istringstream para_istrm(inchar);
		para_istrm >> input;
	}
	else{
		throw std::invalid_argument("Invalid argument type. int");
		//std::cout<<"Error "<<std::endl;//Error message
	}	
}

void read_input_to_size_t(char inchar[], size_t &input){
	if (isdigit(inchar[0])){
		std::istringstream para_istrm(inchar);
		para_istrm >> input;
	}
	else{
		throw std::invalid_argument("Invalid argument type. size_t");
//		std::cout<<"Error"<<std::endl; //Error message
	}	
}


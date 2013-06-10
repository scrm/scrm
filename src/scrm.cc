#include <iostream>
#include <ctime>
#include "seg.h"
#include "forest.h"
#include "random/random_generator.h"
#include "random/mersenne_twister.h"
#include "random/constant_generator.h"
using namespace scrm; 

#ifndef UNITTEST
int main(int argc, char *argv[]){
	//param user_para;
	//check_and_remove("scrm.log");
	if (argc==1 ){
		//scrm_help help;
		//help.print_help();
		print_help();
	}	//else, proceed

    try {
      time_t start_time = time(0);
	    param user_para(argc, argv);

      //Model model = Model(5);
      Model *model = new Model(user_para);

      //Different random generators.
      MersenneTwister *rg = new MersenneTwister(user_para.random_seed);
      //ConstantGenerator *rg = new ConstantGenerator();
      //FakeRandomGenerator *rg = new FakeRandomGenerator();

      // NORMAL RUN
      std::ofstream tree_file;
      for (size_t rep_i=0;rep_i<user_para.nreps;rep_i++){
	      Forest * forest = new Forest(model, rg);
	      forest->buildInitialTree();
	      
	      if (user_para.total_mut==0){
			
		    tree_file.open (user_para.treefile.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
			tree_file << forest->local_root()->tree_topo_bl <<"\n";  
			tree_file.close();	
		  }

	      if (user_para.seg_bool){
			int total_mut;
			double total_bl= forest->local_root()->length_below();

			if (user_para.total_mut!=0){
				total_mut = user_para.total_mut;
				//give the probability ...
				}
			else{
				total_mut=poisson_rand_var(user_para.theta*total_bl);
			}
			  
			forest->seg_data(user_para.treefile, total_mut);
		  }
			if (user_para.rho > 0.0){
				std::ofstream tree_file;
				tree_file.open (user_para.treefile.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
		      //for (size_t i=2; i <= user_para.ith_change; ++i) {
		        //if ( i % 100 == 0 ) std::cout << "active recs done " << i << std::endl;
				  double nextbase=0;
			  while (nextbase < user_para.nsites ) {

		        nextbase=forest->current_base() +forest->random_generator()->sampleExpo(forest->local_tree_length() * forest->writable_model()->recombination_rate() / 4.0 / forest->writable_model()->population_size());
		        forest->set_current_base(nextbase);    
		        forest->sampleNextGenealogy();
		        tree_file << "["<<floor(nextbase) <<"] "<< forest->local_root()->tree_topo_bl <<"\n";
		      }
		      //tree_file  <<"//\n";
			  tree_file.close();	
			}
			
		    tree_file.open (user_para.treefile.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
			tree_file  <<"//\n";
			tree_file.close();	
			//std::cout << forest->getNodes()->size() << std::endl;
			delete forest;
		}
		
      time_t end_time = time(0);
      
      std::cout << "Simulation took about " << end_time - start_time 
                << " second(s)" << std::endl;
                
      if (user_para.log_bool){          
				std::ofstream log_file;
				log_file.open (user_para.log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
				log_file << "Simulation took about " << end_time - start_time << " second(s) \n";
				log_file << "Trees are saved in: "<<user_para.treefile<<"\n";
				log_file.close();
				string system_cmd="cat "+user_para.log_NAME;
				system(system_cmd.c_str());
			}
    }
    catch (const exception &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}
#endif

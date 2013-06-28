#include <iostream>
#include <ctime>
#include "forest.h"
#include "seg.h"
#include "random/random_generator.h"
#include "random/mersenne_twister.h"
#include "random/constant_generator.h"
using namespace scrm; 

#ifndef UNITTEST
int main(int argc, char *argv[]){
  //param user_para;
  //check_and_remove("scrm.log");
  if (argc==1 ){
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
    //
    for (size_t rep_i=0; rep_i<user_para.nreps; ++rep_i){
      Forest * forest = new Forest(model, rg);
      forest->buildInitialTree();
      if (user_para.total_mut==0){	
        std::ofstream tree_file;
        tree_file.open (user_para.treefile.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
        tree_file<<writeTree(forest->local_root(),forest->writable_model()->population_size(), 0.0)<<"\n";
        forest->printTree();
        //std::cout<<writeTree(forest->local_root(),forest->writable_model()->population_size(), 0.0)<<std::endl;
        //std::cout<<writeTree_new(forest->local_root(),forest->writable_model()->population_size())<<std::endl;
        tree_file.close();	
      }

      if (user_para.tmrca_bool){
        std::ofstream tmrca_file;
        tmrca_file.open (user_para.tmrca_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
        tmrca_file << forest->local_root()->height() <<";\n";  
        tmrca_file.close();	
      }
      //vector <seg_data*> seg_data_array;
      seg_data_container * seg_data_array = new seg_data_container(user_para);
      
      int previous_last_for = 1 + min(ceil(forest->next_base()),user_para.nsites)-ceil(forest->current_base());
      if (user_para.rho > 0.0){
        std::ofstream tree_file;
        tree_file.open (user_para.treefile.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
        //string previous_genealogy=writeTree(forest->local_root(),forest->writable_model()->population_size(), 0.0);
         string previous_genealogy=writeTree_new(forest->local_root(),forest->writable_model()->population_size());
        
        while (forest->next_base()< user_para.nsites ) {
          forest->sampleNextGenealogy();
		  seg_data_array->append_new_seg_data(forest);

          string current_genealogy=writeTree_new(forest->local_root(),forest->writable_model()->population_size());
          if (current_genealogy == previous_genealogy){
            previous_last_for=previous_last_for+min(ceil(forest->next_base()),user_para.nsites)-ceil(forest->current_base());
          }
          else{
            tree_file << "["<< previous_last_for <<"] "<< previous_genealogy <<";\n";
            previous_genealogy=current_genealogy;
            previous_last_for=min(ceil(forest->next_base()),user_para.nsites)-ceil(forest->current_base());
          }
        }
		//simulate seg sites ....
        tree_file << "["<< (previous_last_for)  <<"] "<< previous_genealogy <<";\n";
        //tree_file  <<"//\n";
        tree_file.close();	
      }
		seg_data_array->append_new_seg_data(forest);
		seg_data_array->print_to_file();
		delete seg_data_array;

      std::ofstream tree_file;
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

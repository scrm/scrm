#include <iostream>
#include <ctime>

#include "forest.h"
#include "random/random_generator.h"
#include "random/mersenne_twister.h"
#include "random/constant_generator.h"

#ifndef UNITTEST
int main(int argc, char *argv[]){
	//param user_para;
	//check_and_remove("scrm.log");

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
      Forest forest = Forest(model, rg);
      forest.buildInitialTree();

      for (size_t i=2; i <= user_para.ith_change; ++i) {
        if ( i % 100 == 0 ) std::cout << "active recs done " << i << std::endl;
        forest.set_current_base(i);    
        forest.sampleNextGenealogy();
      }

      std::cout << forest.getNodes()->size() << std::endl;

      time_t end_time = time(0);
      
      std::cout << "Simulation took about " << end_time - start_time 
                << " second(s)" << std::endl;
                
      if (user_para.log_bool){          
				std::ofstream log_file;
				log_file.open (user_para.log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
				log_file << "Simulation took about " << end_time - start_time << " second(s) \n";
				log_file.close();
			}
    }
    catch (const exception &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}
#endif

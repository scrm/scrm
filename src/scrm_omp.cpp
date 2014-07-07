
#include <iostream>
#include <ctime>

#include "param.h"
#include "forest.h"
#include "random/random_generator.h"
#include "random/mersenne_twister.h"
#include <climits> // INT_MAX


#ifndef UNITTEST
int main(int argc, char *argv[]){
  try {
    // Organize output
    std::ostream *output = &std::cout;

    Param user_para(argc, argv);

    Model top_model;
    user_para.parse(top_model);

    // Print help if user asked for it
    if (user_para.help()) {
      user_para.printHelp(*output); 
      return EXIT_SUCCESS;
    }
    if (user_para.version()) {
      *output << "scrm " << VERSION << std::endl;
      return EXIT_SUCCESS;
    }

    MersenneTwister rg = MersenneTwister(user_para.random_seed());

    //*output << "scrm " << VERSION << " |" << user_para << std::endl;
    *output << user_para << std::endl;
    *output << rg.seed() << std::endl;

    // Loop over the independent samples
    #pragma omp parallel for schedule(dynamic) 
    for (size_t rep_i=0; rep_i < top_model.loci_number(); ++rep_i) {
      size_t new_seed = (size_t) rg.sampleInt( INT_MAX );
      RandomGenerator* new_rg = new MersenneTwister( new_seed , rg.ff() ); 

      Model* tmp_model = new Model(top_model);
      //cout <<" doing " <<rep_i<<"th rep"<<endl;
      // Mark the start of a new independent sample
      //*output << std::endl << "//" << std::endl;
      ParallelStream mytmp_strm;
      mytmp_strm << "\n//"<<" doing " <<rep_i<<"th rep, thread"<< omp_get_thread_num()<<"\n" ;
      // Now set up the ARG, and sample the initial tree
      Forest forest = Forest(tmp_model, new_rg);
      forest.buildInitialTree();
      //forest.printSegmentSumStats_omp(mytmp_strm);

      while (forest.next_base() < top_model.loci_length()) { 
        // Sample next genealogy
        forest.sampleNextGenealogy();
        //forest.printSegmentSumStats_omp(mytmp_strm);
      }

      //forest.printLocusSumStats_omp(mytmp_strm);
      *output << mytmp_strm.toString();
      delete new_rg;
      delete tmp_model;
    }

    // Clean-up and exit
    rg.clearFastFunc();
    return EXIT_SUCCESS;
  }
  catch (const exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    std::cerr << "Try 'scrm --help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
}
#endif

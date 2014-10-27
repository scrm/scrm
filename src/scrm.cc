/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013, 2014 Paul R. Staab, Sha (Joe) Zhu and Gerton Lunter
 * 
 * This file is part of scrm.
 * 
 * scrm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <iostream>
#include <ctime>

#include "param.h"
#include "forest.h"
#include "random/random_generator.h"
#include "random/mersenne_twister.h"


#ifndef UNITTEST
int main(int argc, char *argv[]){
  try {
    //Organize output
    std::ostream *output = &std::cout;

    // Parse command line arguments
    Param user_para(argc, argv);
    Model model;
    user_para.parse(model);

   std::string tre_str ( "((1:1,2:1):3,3:4);" );
    //std::string tre_str ( "((1:1.1,2:1.1):6,(3:3.3,4:3.3):3.8);" );
    //dout << tre_str << std::endl;
    //std::string tre_str ( "(6:6,(3:3,4:4):5);" );

    MersenneTwister rg = MersenneTwister(1);

    Forest forest = Forest(&model, &rg);

    forest.readNewick ( tre_str );
    
    forest.printLocusSumStats(*output);
          while (forest.next_base() < model.loci_length()) { 
        // Sample next genealogy
        forest.sampleNextGenealogy();
      }
    
     



     ////Print help if user asked for it
    //if (user_para.help()) {
      //user_para.printHelp(*output); 
      //return EXIT_SUCCESS;
    //}
    //if (user_para.version()) {


      //*output << "scrm " << VERSION << std::endl;
      //return EXIT_SUCCESS;
    //}

    //MersenneTwister rg = MersenneTwister(user_para.random_seed());
    //*output << user_para << std::endl;
    //*output << rg.seed() << std::endl;

    //// Create the forest
    //Forest forest = Forest(&model, &rg);

    //// Loop over the independent loci/chromosomes
    //for (size_t rep_i=0; rep_i < model.loci_number(); ++rep_i) {

      //// Mark the start of a new independent sample
      //*output << std::endl << "//" << std::endl;

      //// Now set up the ARG, and sample the initial tree
      //forest.buildInitialTree();
//std::cout  << "contemporaries_.size()"<<forest.contemporaries()->size(0) <<std::endl;
////cout << forest.contemporaries()->size(0) <<std::endl;

      //while (forest.next_base() < model.loci_length()) { 
        //// Sample next genealogy
        //forest.sampleNextGenealogy();
      //}
      
      //forest.printLocusSumStats(*output);
      //forest.clear();
    //}

    // Clean-up and exit
    rg.clearFastFunc();
    return EXIT_SUCCESS;
  }
  catch (const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    std::cerr << "Try 'scrm --help' for more information." << std::endl;
    return EXIT_FAILURE;
  }
}
#endif

/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013 Paul R. Staab, Sha (Joe) Zhu and Gerton Lunter
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
#include "forest.h"
#include "seg.h"
#include "random/random_generator.h"
#include "random/mersenne_twister.h"
#include "random/constant_generator.h"

#ifndef UNITTEST
int main(int argc, char *argv[]){
  if (argc < 3 ){
    std::cout << "Too few command line arguments" << std::endl;
    print_help();
    exit(1);
  }

  try {
    time_t start_time = time(0);
    Param user_para(argc, argv);

    Model* model = user_para.parse();

    MersenneTwister *rg = new MersenneTwister(user_para.random_seed);

    // Organize output
    std::ostream *output = &std::cout;
    if (user_para.log_bool) {
      std::ofstream log_file;
      log_file.open(user_para.log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
      output = &log_file; 
    }

    *output << user_para << std::endl;
    *output << rg->seed() << std::endl;

    //*output << *model;

    /*
       std::ofstream tree_file;
       tree_file.open (user_para.treefile.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
       tree_file  <<"//\n";
       tree_file.close();	
       */

    // Loop over the independent samples
    for (size_t rep_i=0; rep_i < model->loci_number(); ++rep_i) {

      // Mark the start of a new independent sample
      *output << std::endl << "//" << std::endl;

      // Set up a buffer to hold the tree representations
      std::ostringstream tree_buffer;

      // Now set up the ARG, and sample the initial tree
      Forest *forest = new Forest(model, rg);
      forest->buildInitialTree();

      // Set up a buffer to hold the data for segregating sites
      SegDataContainer *seg_data_array = new SegDataContainer(&user_para, forest);

      // Optionally output the TMRCA of the initial coalescent tree in a file
      if (user_para.tmrca_bool()){
        std::ofstream tmrca_file;
        tmrca_file.open (user_para.tmrca_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
        tmrca_file << forest->local_root()->height() <<"\n";  
        tmrca_file.close();	
      }

      // Just output a single tree if the recombination rate is 0
      if (forest->model().mutation_exact_number() == -1 && model->recombination_rate() == 0.0){	
        tree_buffer << writeTree_new(forest->local_root(), forest->model().population_size()) << ";\n";
        seg_data_array->append_new_seg_data(forest);
      }

      // Start main loop, if the recombination rate is nonzero
      if (model->recombination_rate() > 0.0){

        size_t previous_recombination_event;
        size_t next_recombination_event;
        size_t next_event;
        size_t distance_between_events;

        do {

          previous_recombination_event = ceil(forest->current_base());
          next_recombination_event = ceil(forest->next_base());
          next_event = min(1+forest->model().loci_length(), next_recombination_event);    // add 1 since we start at base 1.0
          distance_between_events = next_event - previous_recombination_event;

          // Obtain string representation of current tree
          string previous_genealogy;
          if (user_para.tree_bool) { 
            previous_genealogy = writeTree_new(forest->local_root(),forest->writable_model()->population_size());
          }

          // Sample next genealogy
          forest->sampleNextGenealogy();

          // Sample and store segregating sites data
          seg_data_array->append_new_seg_data(forest);

          // Store current local tree and distance between recombinations in tree buffer
          if (user_para.tree_bool) {
            tree_buffer << "[" << distance_between_events << "] "<< previous_genealogy << ";\n";
          }

        } while ( next_recombination_event == next_event ); // stop if the next recombination event is outside the locus

      }

      if (user_para.tree_bool) {
        *output << tree_buffer.str();
      }

      *output << *seg_data_array;

      delete seg_data_array;
      delete forest;
    }

    //tree_file.close();	
    time_t end_time = time(0);

    std::cout << "Simulation took about " << end_time - start_time 
        << " second(s)" << std::endl;

    if (user_para.log_bool){          
      std::ofstream log_file;
      log_file.open (user_para.log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
      log_file << "Simulation took about " << end_time - start_time << " second(s) \n";
      log_file << "Trees are saved in: " << user_para.tree_NAME << "\n";
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

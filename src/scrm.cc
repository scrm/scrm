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
    Param user_para(argc, argv);

    Model model;
    user_para.parse(model);

    MersenneTwister rg = MersenneTwister(user_para.random_seed());

    // Organize output
    std::ostream *output = &std::cout;

    *output << user_para << std::endl;
    *output << rg.seed() << std::endl;

    double tmrca, tot_bl;

    // Set up a buffer to hold the tree representations
    std::ostringstream tree_buffer;

    // Loop over the independent samples
    for (size_t rep_i=0; rep_i < model.loci_number(); ++rep_i) {

      // Mark the start of a new independent sample
      *output << std::endl << "//" << std::endl;

      // Now set up the ARG, and sample the initial tree
      Forest forest = Forest(&model, &rg);
      forest.buildInitialTree();

      // Set up a buffer to hold the data for segregating sites
      SegDataContainer seg_data_array = SegDataContainer(&user_para, &forest);

      // Just output a single tree if the recombination rate is 0
      if (model.mutation_exact_number() == -1 && model.recombination_rate() == 0.0){	
        tree_buffer << forest.writeTree(forest.local_root()) << ";\n";
        //tree_buffer << writeTree_new(forest.local_root(), forest.model().default_pop_size) << ";\n";
        seg_data_array.append_new_seg_data(&forest);
        tmrca = forest.tmrca();
        tot_bl = forest.tot();
      }

      int i = 0;
      // Start main loop, if the recombination rate is nonzero
      if (model.recombination_rate() > 0.0){

      while (forest.next_base() < model.loci_length()) { 
          // Obtain string representation of current tree
          string previous_genealogy;
          if (user_para.tree_bool) { 
            previous_genealogy = forest.writeTree(forest.local_root());
          }
          tmrca = forest.tmrca();
          tot_bl = forest.tot();

          // Sample next genealogy
          forest.sampleNextGenealogy();

          // Sample and store segregating sites data
          seg_data_array.append_new_seg_data(&forest);

          // Store current local tree and distance between recombinations in tree buffer
          if (user_para.tree_bool && forest.calcSegmentLength() > 0) {
            tree_buffer << "[" << forest.calcSegmentLength() << "]" << previous_genealogy << ";\n";
          }
        }

      }

      if (user_para.tree_bool) {
        *output << tree_buffer.str();
      }

      if (user_para.tmrca_bool){
        *output << "time:\t"<<tmrca<< "\t"<<tot_bl <<"\n";  
      }

      *output << seg_data_array;
    }

  }
  catch (const exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}
#endif

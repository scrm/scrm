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

  //try {
    time_t start_time = time(0);
    Param user_para(argc, argv);

    Model* model = user_para.parse();

    //Different random generators.
    MersenneTwister *rg = new MersenneTwister(user_para.random_seed);
    //ConstantGenerator *rg = new ConstantGenerator();
    //FakeRandomGenerator *rg = new FakeRandomGenerator();

    // Organize output
    std::ostream *output = &std::cout;
    if (user_para.log_bool) {
      std::ofstream log_file;
      log_file.open(user_para.log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
      output = &log_file; 
    }

    *output << user_para << std::endl;
    *output << user_para.random_seed << std::endl;
    //*output << *model;

    /*
       std::ofstream tree_file;
       tree_file.open (user_para.treefile.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
       tree_file  <<"//\n";
       tree_file.close();	
       */

    string previous_genealogy;

    for (size_t rep_i=0; rep_i < model->loci_number(); ++rep_i) {
      std::ostringstream tree_buffer;
      *output << std::endl << "//" << std::endl;

      Forest *forest = new Forest(model, rg);
      forest->buildInitialTree();

      if (user_para.tree_bool && forest->model().mutation_exact_number() == -1 && model->recombination_rate() == 0.0){	
        tree_buffer << writeTree_new(forest->local_root(), forest->model().population_size()) << ";\n";
      }

      if (user_para.tmrca_bool()){
        std::ofstream tmrca_file;
        tmrca_file.open (user_para.tmrca_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
        tmrca_file << forest->local_root()->height() <<"\n";  
        tmrca_file.close();	
      }

      SegDataContainer *seg_data_array = new SegDataContainer(&user_para, forest);

      int previous_last_for = 1 + min((size_t)ceil(forest->next_base()), forest->model().loci_length())-ceil(forest->current_base());
      if (model->recombination_rate() > 0.0){
        if (user_para.tree_bool) { 
          previous_genealogy=writeTree_new(forest->local_root(),forest->writable_model()->population_size());
        }

        while ( forest->next_base() < model->loci_length() ) {
          forest->sampleNextGenealogy();
          seg_data_array->append_new_seg_data(forest);

          if (user_para.tree_bool) {
            tree_buffer << "["<< previous_last_for <<"] "<< previous_genealogy <<";\n";
            previous_last_for=min((size_t)ceil(forest->next_base()), forest->model().loci_length())-ceil(forest->current_base());
          }
        }		
        if (user_para.tree_bool) tree_buffer << "["<< (previous_last_for)  <<"] "<< previous_genealogy <<";\n";

      }
      seg_data_array->append_new_seg_data(forest);

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
      //string system_cmd="cat "+user_para.log_NAME;
      //system(system_cmd.c_str()); //VERY UGLY -Paul
    }
  //}
  /*
  catch (const exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  */
}
#endif

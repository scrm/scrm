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

#include "orientedForest.h"

void OrientedForest::calculate(const Forest &forest) {
  if (forest.model().recombination_rate() == 0.0) {
    output_buffer_ << generateTree(forest) << ";\n";  
  } else {
    if (forest.calcSegmentLength(forest.model().finite_sites()) == 0.0) return;
    output_buffer_ << "[" << forest.calcSegmentLength(forest.model().finite_sites()) << "]" 
                   << generateTree(forest) << ";\n";  
  }
}

void OrientedForest::printLocusOutput(std::ostream &output) {
  output << output_buffer_.str();  
  output_buffer_.str("");
  output_buffer_.clear();
}

/**
 * @brief Prints a part of the tree in newick format
 *
 * @param node The root of the subtree that will be printed
 *
 * @return A part of the tree in newick format
 */
std::string OrientedForest::generateTree(const Forest &forest) {
  this->init_OF_label( forest );

  // a node can either be in sample (labelled), or not in sample, 
  // if it is insample, check if its parent is labelled.
  for (auto it = forest.getNodes()->iterator(); it.good(); ++it) {
    if ( !(*it)->local() ) continue;
    Node* tmp_parent = (*it)->parent();
    if ( (*it)->in_sample() && (*it)->parent_OF_label() == 0 ) {
        if ( !tmp_parent->is_root() ) this->tmp_label++;
        //tmp_parent->set_OF_label( tmp_parent->is_root() ? 0:this->tmp_label );
        tmp_parent->set_OF_label( tmp_parent->is_root() ? 0:this->tmp_label );
        cout << (*it)->OF_label() << "("<< (*it)->parent_OF_label() << ")" << ",";
    }
    else if (  (*it)->in_sample()  ){
        cout << (*it)->OF_label() << "("<< (*it)->parent_OF_label() << ")" << ",";
        }
    else if ( (*it)->parent_OF_label() == 0 ){ // (*it)->in_sample() 
        if ( !tmp_parent->is_root() ) this->tmp_label++;
        tmp_parent->set_OF_label( tmp_parent->is_root() ? 0:this->tmp_label );
        cout << (*it)->OF_label() << "("<< (*it)->parent_OF_label() << ")" << ",";
        }
    else {
        cout << (*it)->OF_label() << "("<< (*it)->parent_OF_label() << ")" << ",";
        }
    }
    cout<<endl;

  return std::string("ok");
}


void OrientedForest::init_OF_label( const Forest &forest ){
  this->OF_labels.clear();
  //forest.reset_OF_label();
  for (auto it = forest.getNodes()->iterator(); it.good(); ++it) {
    if ( !(*it)->local() ) continue;
    //(*it)->set_OF_label( (*it)->in_sample() ? (*it)->label():0 );
    //(*it)->set_parent_OF_label( (size_t)0 );
    this->OF_labels.push_back( (size_t)0);    
  }
    // Initialize the first parent node to tip node 1, as the number of tips plus one.
  this->tmp_label = forest.model().sample_size();
}

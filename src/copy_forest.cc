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

//copy_forest.cc
#include "forest.h"
Forest::Forest(
Forest * current_forest /*! Forest that needs to be duplicated */, 
bool entire_ARG /*! if entire_ARG is false, only copy the local tree (McVean2005) */
){ 
	//Initialize 
	this->nodes_ = new NodeContainer();
	//this->set_local_root(NULL);
	//this->set_primary_root(NULL);
  
  for(size_t i = 0; i < this->getNodes()->size(); ++i) {
	  this->nodes()->get_copy(i)->index = i;  
  }

	dout<<"start adding nodes"<<std::endl;
	dout<<current_forest->getNodes()->size()<<endl;
	size_t reverse_i=0;
	for(; reverse_i < current_forest->getNodes()->size(); reverse_i++) {
		if (current_forest->getNodes()->get(reverse_i)->label()==0){
		  break;
		}
	}
	dout<<reverse_i<<" tips"<<endl;
	

	for(size_t i = reverse_i; i > 0; i--) {
	  double height = current_forest->getNodes()->get(i-1)->height();
	  bool local = current_forest->getNodes()->get(i-1)->local();
	  size_t last_update = current_forest->getNodes()->get(i-1)->last_update();
	  size_t samples_below = current_forest->getNodes()->get(i-1)->samples_below();
	  double length_below = current_forest->getNodes()->get(i-1)->length_below();
	  int label = current_forest->getNodes()->get(i-1)->label();
	  Node * new_node=new Node(height,local,last_update, samples_below, length_below, label);
	  new_node->index = current_forest->getNodes()->get(i-1)->index;
	  this->nodes()->add(new_node);
	  
	}
	for(size_t i=reverse_i; i < current_forest->getNodes()->size(); i++) {
	  double height = current_forest->getNodes()->get(i)->height();
	  bool local = current_forest->getNodes()->get(i)->local();
	  size_t last_update = current_forest->getNodes()->get(i)->last_update();
	  size_t samples_below = current_forest->getNodes()->get(i)->samples_below();
	  double length_below = current_forest->getNodes()->get(i)->length_below();
	  int label = current_forest->getNodes()->get(i)->label();
	  Node * new_node=new Node(height,local,last_update, samples_below, length_below, label);
	  new_node->index = current_forest->getNodes()->get(i)->index;
	  this->nodes()->add(new_node);	  
	}	
	
	for(size_t i = 0; i < current_forest->getNodes()->size(); i++) {
		if (current_forest->nodes()->get_copy(i)->parent()!=NULL){
			this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_parent(this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->parent()->index));		
		}
		else{
			this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_parent(NULL);
		}
		
		if (current_forest->nodes()->get_copy(i)->first_child()!=NULL){
			this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_first_child(this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->first_child()->index));
		}
		else{
			this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_first_child(NULL);
		}
		
		if (current_forest->nodes()->get_copy(i)->second_child()==NULL){
			this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_second_child(NULL);
			
		}
		else{
			this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_second_child(this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->second_child()->index));
		}
	}
	this->set_local_root(this->nodes()->get_copy(current_forest->local_root()->index));
	this->set_primary_root(this->nodes()->get_copy(current_forest->primary_root()->index));
	this->set_model(current_forest->writable_model());
	this->set_sample_size(current_forest->sample_size());
	this->set_current_base(current_forest->current_base());
	this->set_random_generator(current_forest->random_generator());
	this->set_expo_sample(current_forest->expo_sample());
	this->set_prune_countdown(current_forest->prune_countdown());
	this->set_pruning(current_forest->pruning());
	this->set_next_base();

	dout<<"  #################### check copied forest ###############"<<std::endl;
	assert(this->printTree());
	assert(this->checkTree());
	dout<<"  #################### check copied forest finished ###############"<<std::endl<<std::endl;
}


//copy_forest.cc
#include "forest.h"
Forest::Forest(Forest * current_forest, bool entire_ARG){ // if entire_ARG is false, only copy the local tree (McVean2005)
	//this->initialize(current_forest->writable_model() ,current_forest->random_generator());
	//current_forest->printTree();
	//dout<<"There are "<<current_forest->getNodes()->size()<<" nodes"<<endl;
	this->nodes_ = new NodeContainer();
	dout<<"start adding nodes"<<std::endl;
	
	for(size_t i = 0; i < current_forest->getNodes()->size(); i++) {
	  Node * new_node=new Node(current_forest->getNodes()->get(i)->height(), current_forest->getNodes()->get(i)->local(), current_forest->getNodes()->get(i)->last_update(), current_forest->getNodes()->get(i)->samples_below(), current_forest->getNodes()->get(i)->length_below(),current_forest->getNodes()->get(i)->label());
	  this->nodes()->add(new_node);
	}
	//dout<<"finished adding nodes "<< this->getNodes()->size()<<std::endl;

	for(size_t i = 0; i < this->getNodes()->size(); i++) {
		//dout<<"i = "<<i<<endl;
		if (current_forest->nodes()->get_copy(i)->parent()!=NULL){
			//cout<<"it is not null"<<endl;
			this->nodes()->get_copy(i)->set_parent(this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->parent()->index));		
		}
		else{
			//cout<<"it is null"<<endl;
			this->nodes()->get_copy(i)->set_parent(NULL);
		}
	  //dout<<"added parent"<<endl;
		
		if (current_forest->nodes()->get_copy(i)->first_child()!=NULL){
			//cout<<"it is not null"<<endl;

			this->nodes()->get_copy(i)->set_first_child(this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->first_child()->index));
		}
		else{
//cout<<"it is null"<<endl;
			this->nodes()->get_copy(i)->set_first_child(NULL);
		}
	  //dout<<"added first_child"<<endl;
		
		if (current_forest->nodes()->get_copy(i)->second_child()!=NULL){
			//cout<<"it is not null"<<endl;
			this->nodes()->get_copy(i)->set_second_child(this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->second_child()->index));
		}
		else{
//cout<<"it is null"<<endl;
			this->nodes()->get_copy(i)->set_second_child(NULL);
		}
	  //dout<<"added second_child"<<endl;

	}

	dout<<"finished linking nodes"<<std::endl;
	
	this->set_model(current_forest->writable_model());
	this->set_sample_size(current_forest->sample_size());
	this->set_current_base(current_forest->current_base());
  this->set_random_generator(current_forest->random_generator());
  this->set_local_root(this->nodes()->get_copy(current_forest->local_root()->index));
  this->set_primary_root(this->nodes()->get_copy(current_forest->primary_root()->index));
  this->set_expo_sample(current_forest->expo_sample());
  this->set_prune_countdown(current_forest->prune_countdown());
  this->set_pruning(current_forest->pruning());

	dout<<"  #################### check copied forest ###############"<<std::endl;
	assert(this->printTree());
  //assert(this->printNodes());
  assert(this->checkTree());
  dout<<"  #################### check copied forest finished ###############"<<std::endl<<std::endl;
}


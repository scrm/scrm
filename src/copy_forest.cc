//copy_forest.cc
#include "forest.h"
Forest::Forest(Forest * current_forest, bool entire_ARG){ // if entire_ARG is false, only copy the local tree (McVean2005)
	//this->initialize(current_forest->writable_model() ,current_forest->random_generator());
	//current_forest->printTree();
	//dout<<"There are "<<current_forest->getNodes()->size()<<" nodes"<<endl;
	
	//Initialize 
	this->nodes_ = new NodeContainer();
	//this->nodes_ = current_forest->nodes();

	//this->set_local_root(NULL);
	//this->set_primary_root(NULL);
	dout<<"start adding nodes"<<std::endl;
	
	dout<<current_forest->getNodes()->size()<<endl;
	
	//Node * previous=NULL;
	size_t reverse_i=0;
	for(; reverse_i < current_forest->getNodes()->size(); reverse_i++) {
		if (current_forest->getNodes()->get(reverse_i)->label()==0){
		  break;
		}
	}
	cout<<reverse_i<<" tips"<<endl;
	

	for(size_t i = reverse_i; i > 0; i--) {
	  //cout<<i<<endl;
	 	 double height = current_forest->getNodes()->get(i-1)->height();
	  bool local = current_forest->getNodes()->get(i-1)->local();
	  size_t last_update = current_forest->getNodes()->get(i-1)->last_update();
	  size_t samples_below = current_forest->getNodes()->get(i-1)->samples_below();
	  double length_below = current_forest->getNodes()->get(i-1)->length_below();
	  int label = current_forest->getNodes()->get(i-1)->label();
	  //Node(double height, bool local, size_t last_update, size_t samples_below, double length_below, int label);
	  Node * new_node=new Node(height,local,last_update, samples_below, length_below, label);
	  new_node->index = current_forest->getNodes()->get(i-1)->index;
	  this->nodes()->add(new_node);
	  
	}
	for(size_t i=reverse_i; i < current_forest->getNodes()->size(); i++) {
	//cout<<i<<endl;
	  double height = current_forest->getNodes()->get(i)->height();
	  bool local = current_forest->getNodes()->get(i)->local();
	  size_t last_update = current_forest->getNodes()->get(i)->last_update();
	  size_t samples_below = current_forest->getNodes()->get(i)->samples_below();
	  double length_below = current_forest->getNodes()->get(i)->length_below();
	  int label = current_forest->getNodes()->get(i)->label();
	  //Node(double height, bool local, size_t last_update, size_t samples_below, double length_below, int label);
	  Node * new_node=new Node(height,local,last_update, samples_below, length_below, label);
	  new_node->index = current_forest->getNodes()->get(i)->index;
	  this->nodes()->add(new_node);	  
	}	
	
	




	for(size_t i = 0; i < current_forest->getNodes()->size(); i++) {
		//dout<<i<<"th node"<<endl;
		if (current_forest->nodes()->get_copy(i)->parent()!=NULL){
			//cout<<"parent is not null"<<endl;
			this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_parent(this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->parent()->index));		
		}
		else{
			//cout<<"parent is null"<<endl;
			this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_parent(NULL);
		}
	  //dout<<"added parent"<<endl;
		
		if (current_forest->nodes()->get_copy(i)->first_child()!=NULL){
			//cout<<"it is not null"<<endl;

			this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_first_child(this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->first_child()->index));
		}
		else{
//cout<<"it is null"<<endl;
			this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_first_child(NULL);
		}
	  //dout<<"added first_child"<<endl;
		
		if (current_forest->nodes()->get_copy(i)->second_child()==NULL){
			//cout<<"it is null"<<endl;
			this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_second_child(NULL);
			
		}
		else{
//cout<<"it is not null"<<endl;
			this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_second_child(this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->second_child()->index));
		}
	  //dout<<"added second_child"<<endl;
	  
		//if (i== (current_forest->getNodes()->size()-1)){
			//this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_next(NULL);
			
		//}
		//else{
			//this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_next(this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->next()->index));
		//}
		////assert(current_forest->nodes()->get_copy(i)->previous());
		//if (i==0){
			//this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_previous(NULL);

		//}
		//else{
			////cout<<"something wrong here"<<endl;
			//this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->index)->set_previous(this->nodes()->get_copy(current_forest->nodes()->get_copy(i)->previous()->index));
			
		//}
	  //dout<<"added second_child"<<endl;

	}
	this->set_local_root(this->nodes()->get_copy(current_forest->local_root()->index));
	this->set_primary_root(this->nodes()->get_copy(current_forest->primary_root()->index));
	
	
	//assert(current_forest->printNodes());
	//assert(this->printNodes());
	dout<<"finished adding nodes "<< this->getNodes()->size()<<std::endl;
	cout<<"i  old index    new index"<<endl;
	for(size_t i = 0; i < current_forest->getNodes()->size(); i++) {
		cout<<i<<"  "<<current_forest->nodes()->get_copy(i)->index<<" "<<current_forest->nodes()->get_copy(i)->parent_height();
		cout<<"   "<<this->nodes()->get_copy(i)->index<<"  "<<this->nodes()->get_copy(i)->parent_height();
		cout<<endl;
	}
	
	
	dout<<"finished linking nodes"<<std::endl;
		this->set_model(current_forest->writable_model());
	this->set_sample_size(current_forest->sample_size());
	this->set_current_base(current_forest->current_base());
	this->set_random_generator(current_forest->random_generator());
	this->set_expo_sample(current_forest->expo_sample());
	this->set_prune_countdown(current_forest->prune_countdown());
	this->set_pruning(current_forest->pruning());


	dout<<"  #################### check copied forest ###############"<<std::endl;
	assert(this->printTree());
	//assert(this->printNodes());
	assert(this->checkTree());
	dout<<"  #################### check copied forest finished ###############"<<std::endl<<std::endl;
}


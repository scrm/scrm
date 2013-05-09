//copy_forest.cc
#include"copy_forest.h"
//NodeContainer

Forest duplicate_forest(Forest * current_forest){

	Forest new_forest;
  //new_forest.nodes_=current_forest->getNodes();
	for(size_t i = 0; i < current_forest->getNodes()->size(); ++i) {
		Node * new_node=new Node;
		// copy everything for newnode from current_forest->getNodes()->get(i)
	//new_node->  = current_forest->getNodes()->get(i)->
  new_node->set_height(current_forest->getNodes()->get(i)->height);
  new_node->set_local(current_forest->getNodes()->get(i)->local);
  //this->set_parent(current_forest->getNodes()->get(i)->NULL);
  //this->set_second_child(current_forest->getNodes()->get(i)->NULL);
  //this->set_first_child(current_forest->getNodes()->get(i)->NULL);
  new_node->set_last_update(current_forest->getNodes()->get(i)->last_update);
  new_node->set_samples_below(current_forest->getNodes()->get(i)->samples_below);
  new_node->set_length_below(current_forest->getNodes()->get(i)->length_below);
  //this->set_previous(current_forest->getNodes()->get(i)->NULL);
  //this->set_next(current_forest->getNodes()->get(i)->NULL);	
	}

	new_forest.set_model(current_forest->model());
  
  
  
  new_forest.set_sample_size(current_forest->sample_size());
  
  new_forest.set_current_base(current_forest->current_base());
  
  //this->nodes_ = new NodeContainer();
  
  new_forest.set_random_generator(current_forest->rg);
  
  
  
  //new_forest.set_local_root(current_forest->local_root());
  
  //new_forest.set_primary_root(current_forest->primary_root());
  
  //expo_sample_;      // Placeholder for exp(1) sampled values
  //size_t prune_countdown_;  // We will prune once this countdown reaches 0
  //bool pruning_;

  //RandomGenerator* random_generator_;

  


  //NodeContainer *nodes() { return this->nodes_; }

	
	
	return new_forest;
}


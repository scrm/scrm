//copy_forest.cc
#include"copy_forest.h"

Forest duplicate_forest(Forest  current_forest){
	Forest new_forest;
	
	//dout<<"start adding nodes"<<std::endl;
	
	for(size_t i = 0; i < current_forest.getNodes()->size(); ++i) {
	  Node * new_node=new Node;
		// copy everything for newnode from current_forest.getNodes()->get(i)
		
		new_node->set_height(current_forest.getNodes()->get(i)->height());
	  
	  new_node->set_local(current_forest.getNodes()->get(i)->local());
	  
	  new_node->set_last_update(current_forest.getNodes()->get(i)->last_update());
	  new_node->set_samples_below(current_forest.getNodes()->get(i)->samples_below());
	  new_node->set_length_below(current_forest.getNodes()->get(i)->length_below());
	  new_forest.nodes()->add(new_node);
	}
	//dout<<"finished adding nodes"<<std::endl;

	for(size_t i = 0; i < new_forest.getNodes()->size(); ++i) {
		//dout<<"i = "<<i<<endl;
		if (current_forest.nodes()->get_copy(i)->parent()!=NULL){
			new_forest.nodes()->get_copy(i)->set_parent(new_forest.nodes()->get_copy(current_forest.nodes()->get_copy(i)->parent()->index));		
		}
		else{
			new_forest.nodes()->get_copy(i)->set_parent(NULL);
		}
	  //dout<<"added parent"<<endl;
		
		if (current_forest.nodes()->get_copy(i)->first_child()!=NULL){
			new_forest.nodes()->get_copy(i)->set_first_child(new_forest.nodes()->get_copy(current_forest.nodes()->get_copy(i)->first_child()->index));
		}
		else{
			new_forest.nodes()->get_copy(i)->set_first_child(NULL);
		}
	  //dout<<"added first_child"<<endl;
		
		if (current_forest.nodes()->get_copy(i)->second_child()==NULL){
			new_forest.nodes()->get_copy(i)->set_second_child(NULL);
		}
		else{
			new_forest.nodes()->get_copy(i)->set_second_child(new_forest.nodes()->get_copy(current_forest.nodes()->get_copy(i)->second_child()->index));
		}
	  //dout<<"added second_child"<<endl;

	}

	//dout<<"finished linking nodes"<<std::endl;
	
	new_forest.set_model(current_forest.writable_model());
	new_forest.set_sample_size(current_forest.sample_size());
	new_forest.set_current_base(current_forest.current_base());
  new_forest.set_random_generator(current_forest.random_generator());
  new_forest.set_local_root(new_forest.nodes()->get_copy(current_forest.local_root()->index));
  new_forest.set_primary_root(new_forest.nodes()->get_copy(current_forest.primary_root()->index));
  new_forest.set_expo_sample(current_forest.expo_sample());
  new_forest.set_prune_countdown(current_forest.prune_countdown());
  new_forest.set_pruning(current_forest.pruning());

	dout<<"  #################### check copied forest ###############"<<std::endl;
	assert(new_forest.printTree());
  assert(new_forest.printNodes());
  assert(new_forest.checkTree());
  dout<<"  #################### check copied forest finished ###############"<<std::endl<<std::endl;
	return new_forest;
}


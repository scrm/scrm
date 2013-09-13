#include"printtree.h"


std::string writeTree(Node * node, Model * model){
	if(node->first_child() == NULL && ((node->label())>0)){ // real tip node
		std::ostringstream label_strm;
		label_strm<<node->label();
		return label_strm.str();
	}
	else{
		double *t1 = new double;
		*t1 = 0;
		Node *left = tracking_local_node_with_bl(node->first_child(), t1, model);
		//double t1=node->height()- left->height();
		std::ostringstream t1_strm;
		t1_strm << *t1;
		double *t2 = new double;
		*t2 = 0;
		Node *right = tracking_local_node_with_bl(node->second_child(), t2, model);
		//double t2=node->height()- right->height();
		std::ostringstream t2_strm;
		t2_strm << *t2;
		return "("+writeTree(left,model)+":"+t1_strm.str()+","+ writeTree(right,model)+":"+t2_strm.str() +")";
	}
}


Node * tracking_local_node_with_bl(Node * node, double * bl, Model * model){
  assert( node->local() );
  size_t pop1 = model->getPopLayerIatHeight(node->parent_height());
  size_t pop2 = model->getPopLayerIatHeight(node->height()); 
	//double pop1HeightBottom = model->change_times_[pop1];
	//double pop2HeightTop =  model->change_times_[pop2+1];

  //dout<<"*****"<<std::endl;
  //dout << node->parent() << " is in Population "<< pop1<<" of size "<< model->popSizeAtLayerI(pop1, node->parent()->population())<<" "<<(node->parent_height() - pop1HeightBottom ) /4/ model->popSizeAtLayerI(pop1,node->parent()->population()) <<std::endl;//	CHECK
  //dout<<std::endl<<"pop2HeightTop="<<pop2HeightTop<<std::endl;
  //dout << node << " is in Population "<< pop2<<" of size "<< model->popSizeAtLayerI(pop2, node->population())<<" ("<<(pop2HeightTop - node->height() )<<") "<<(pop2HeightTop - node->height() ) /4/ model->popSizeAtLayerI(pop2, node->population());//	CHECK
  double dummy = 0;
  
  // condition is not necessary
  if (pop1 == pop2){
	  dummy += (node->parent_height() - node->height() ) /4/ model->popSizeAtLayerI(pop1,node->population());
  }
  else{
	double pop1HeightBottom = model->change_times_[pop1];
	double pop2HeightTop =  model->change_times_[pop2+1];
	for (size_t i = pop2+1; i<pop1; i++ ){
		dummy += (model->change_times_[i+1]- model->change_times_[i]) /4/ model->default_pop_size;
	}

	dummy += (node->parent_height() - pop1HeightBottom ) /4/ model->default_pop_size 
			+ (pop2HeightTop - node->height() ) /4/ model->default_pop_size	;
  }
  
  
  *bl += dummy;
  
  if (node->in_sample()){
    dout<<", is tip node, it is local, return."<<std::endl;
    return node;
  }

  assert( node->first_child() != NULL );
  if ( node->second_child() == NULL ){
    return tracking_local_node_with_bl(node->first_child(), bl, model);
  }

  else if (node->first_child()->local() && node->second_child()->local()){
    dout<< ", is an internal local node, return"<<std::endl;
    return node;
  }
  
  else if (!node->first_child()->local() ){ 
    assert( node->second_child()->local() );
    dout<< ", is internal node with one non-local branch" << std::endl;
    return tracking_local_node_with_bl(node->second_child(), bl, model);
  }

  else{
    assert( node->first_child()->local() );
    dout<< ", is internal node with one non-local branch" << std::endl;
    return tracking_local_node_with_bl(node->first_child(), bl, model);
  }
}


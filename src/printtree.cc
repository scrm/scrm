#include "forest.h"
#include <sstream>

/*! \todo CHECK!!!! NEWICK FORMAT STRING TMRCA AND BL ARE INCORRECT!!!*/

std::string Forest::writeTree(Node * node /*!< Root of the subtree string*/){
	//if(node->first_child() == NULL && ((node->label())>0)){ // real tip node
	if (node->in_sample()){
		std::ostringstream label_strm;
		label_strm<<node->label();
		return label_strm.str();
	}
	else{
		double *t1 = new double;
		*t1 = 0;
		Node *left = this->tracking_local_node_with_bl(node->first_child(), t1);
		std::ostringstream t1_strm;
		t1_strm << *t1;
		double *t2 = new double;
		*t2 = 0;
		Node *right = this->tracking_local_node_with_bl(node->second_child(), t2);
		std::ostringstream t2_strm;
		t2_strm << *t2;
		return "("+this->writeTree(left)+":"+t1_strm.str()+","+ this->writeTree(right)+":"+t2_strm.str() +")";
	}
}

void Forest::tracking_bl(Node * node, double * bl){
  size_t pop1 = model_->getPopLayerIatHeight(node->parent_height());
  size_t pop2 = model_->getPopLayerIatHeight(node->height()); 
  double dummy = 0;
  
  if (pop1 == pop2){
	  dummy += (node->parent_height() - node->height() ) / 4 / model_->popSizeAtLayerI(pop1, node->population());
  }
  else{
	double pop1HeightBottom = model_->change_times_[pop1];
	double pop2HeightTop =  model_->change_times_[pop2+1];
	for (size_t i = pop2 + 1; i < pop1; i++ ){
		dummy += (model_->change_times_[i+1]- model_->change_times_[i]) / 4 / model_->default_pop_size;
	}
	dummy += (node->parent_height() - pop1HeightBottom ) / 4 / model_->default_pop_size 
			+ (pop2HeightTop - node->height() ) / 4 / model_->default_pop_size	;
  }
  
  *bl += dummy;
	
}

Node * Forest::tracking_local_node_with_bl(Node * node, double * bl){
  assert( node->local() );
 
  this->tracking_bl(node,bl);
  
  if (node->in_sample()){
    dout<<", is tip node, it is local, return."<<std::endl;
    return node; /*!< External node */
  }

  assert( node->first_child() != NULL );
  if ( node->second_child() == NULL ){
    return this->tracking_local_node_with_bl(node->first_child(), bl); /*!< Internal node with one child. Keep on searching */
  }

  else if (node->first_child()->local() && node->second_child()->local()){
    dout<< ", is an internal local node, return"<<std::endl;
    return node; /*!< Internal node found */
  }
  
  else if (!node->first_child()->local() ){ 
    assert( node->second_child()->local() );
    dout<< ", is internal node with one non-local branch" << std::endl;
    return this->tracking_local_node_with_bl(node->second_child(), bl); /*!< Internal node, one of the child is non-local. Keep on searching */
  }

  else{
    assert( node->first_child()->local() );
    dout<< ", is internal node with one non-local branch" << std::endl;
    return this->tracking_local_node_with_bl(node->first_child(), bl); /*!< Internal node, one of the child is non-local. Keep on searching */
  }
}

double Forest::tmrca(){/*! \todo Still working progress*/
	double *t1 = new double;
	*t1 = 0;
	//Node *left = this->tracking_local_node_with_bl(local_root_, t1);
	this->tracking_tip_node_with_bl(local_root_->first_child(), t1);
	return *t1;
}


Node * Forest::tracking_tip_node_with_bl(Node * node, double * bl){/*! \todo Still working progress*/
  assert( node->local() );

  this->tracking_bl(node,bl);
  
  if (node->in_sample()){
    dout<<", is tip node, it is local, return."<<std::endl;
    return node; /*!< External node */
  }

  assert( node->first_child() != NULL );
  if ( node->second_child() == NULL ){
    return this->tracking_tip_node_with_bl(node->first_child(), bl); /*!< Internal node with one child. Keep on searching */
  }

  else if (node->first_child()->local() && node->second_child()->local()){
    dout<< ", is an internal local node, return"<<std::endl;
    return this->tracking_tip_node_with_bl(node->first_child(), bl); /*!< Internal node, one of the child is non-local. Keep on searching */
  }
  
  else if (!node->first_child()->local() ){ 
    assert( node->second_child()->local() );
    dout<< ", is internal node with one non-local branch" << std::endl;
    return this->tracking_tip_node_with_bl(node->second_child(), bl); /*!< Internal node, one of the child is non-local. Keep on searching */
  }

  else{
    assert( node->first_child()->local() );
    dout<< ", is internal node with one non-local branch" << std::endl;
    return this->tracking_tip_node_with_bl(node->first_child(), bl); /*!< Internal node, one of the child is non-local. Keep on searching */
  }
}

double Forest::tot_below(Node * node){
	if (node->in_sample()){
		return 0;
	}
	else{
		double *t1 = new double;
		*t1 = 0;
		Node *left = this->tracking_local_node_with_bl(node->first_child(), t1);
		double *t2 = new double;
		*t2 = 0;
		Node *right = this->tracking_local_node_with_bl(node->second_child(), t2);
		return tot_below(left) + *t1 + tot_below(right) + *t2;
	}
}

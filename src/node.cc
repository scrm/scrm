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

#include "node.h"
//#include<iostream>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<iostream>


Node::Node() { init(); };
Node::Node(double height) { init(height, true); };
Node::Node(double height, bool local) { init(height, local); };
Node::Node(double height, bool local, double last_update) { 
  init(height, local, last_update); };
Node::Node(double height, bool local, double last_update, size_t samples_below)
  { init(height, local, last_update, samples_below); };
Node::Node(double height, bool local, double last_update, size_t samples_below, double length_below)
  { init(height, local, last_update, samples_below, length_below); };
Node::Node(double height, bool local, double last_update, size_t samples_below, double length_below, size_t label)
  { init(height, local, last_update, samples_below, length_below, label); };
  
//Node::~Node() {delete[] marginal_likelihood;}
Node::~Node() {};

void Node::init(double height, bool local, double last_update,
                size_t samples_below, double length_below, size_t label) {
  this->set_height(height);
  this->set_local(local);
  this->set_last_update(last_update);
  this->set_samples_below(samples_below);
  this->set_length_below(length_below);
  this->set_label(label);
  this->set_parent(NULL);
  this->set_second_child(NULL);
  this->set_first_child(NULL);
  this->set_previous(NULL);
  this->set_next(NULL);
  this->mutation_state=false;
  this->set_mut_num(0);
  this->set_population(0);
}

Node* Node::parent() const {
  //if (this->parent_ == NULL) 
  //  throw std::logic_error("Trying to access parent of ultimate root"); 
  return this->parent_; 
}

void Node::change_child(Node* from, Node* to) {
  if ( this->first_child() == from )
    this->set_first_child(to);
  else if ( this->second_child() == from )
    this->set_second_child(to);
  else throw std::invalid_argument("Can't find child node to replace");
}

void Node::remove_child(Node* child) {
  if ( this->first_child() == child ) {
    this->set_first_child(this->second_child());
    this->set_second_child(NULL);
    return;
  } 

  if ( this->second_child() == child ) {
    this->set_second_child(NULL);
    return;
  }

  throw std::invalid_argument("Can't find child to delete");
}

bool Node::in_sample() const {
  return ( this->label() != 0 ); 
}

//Node * tracking_local_node(Node * node){ /*! \todo use node->samples_below(), and move it under the forest class*/
  //assert( node->local() );
  //dout << node << std::endl;
  //if (node->in_sample()){
    //dout<<" is tip node, it is local, return."<<std::endl;
    //return node;
  //}

  //assert( node->first_child() != NULL );
  //if ( node->second_child() == NULL ){
    //return tracking_local_node(node->first_child());
  //}

  //else if (node->first_child()->local() && node->second_child()->local()){
    //dout<< " is an internal local node, return"<<std::endl;
    //return node;
  //}
  
  //else if (!node->first_child()->local() ){ 
    //assert( node->second_child()->local() );
    //dout<< "is internal node with one non-local branch" << std::endl;
    //return tracking_local_node(node->second_child());
  //}

  //else{
    //assert( node->first_child()->local() );
    //dout<< "is internal node with one non-local branch" << std::endl;
    //return tracking_local_node(node->first_child());
  //}
//}

//std::string writeTree_new(Node * node, int npop){ /*! \todo TO BE REMOVED */
	//if(node->first_child() == NULL && ((node->label())>0)){ // real tip node
		//std::ostringstream label_strm;
		//label_strm<<node->label();
		//return label_strm.str();
	//}
	//else{
		//Node *left = tracking_local_node(node->first_child());
		//double t1=node->height()- left->height();
		//std::ostringstream t1_strm;
		//t1_strm << t1/4/npop;
		//Node *right = tracking_local_node(node->second_child());
		//double t2=node->height()- right->height();
		//std::ostringstream t2_strm;
		//t2_strm << t2/4/npop;
		//return "("+writeTree_new(left,npop)+":"+t1_strm.str()+","+ writeTree_new(right,npop)+":"+t2_strm.str() +")";
	//}
//}



///**
 //* Extract the string to represent the subtree that is descendant from node.

//*/
//std::string writeTree(Node * node, int npop, double bl_above_parent){ /*! \todo TO BE REMOVED */
	//if (!node->local()){
		//return "("+ writeTree(node->first_child(),npop,0) +","+ writeTree(node->second_child(),npop,0)+")";
		//}
	//else{	
		//if(node->first_child() == NULL && ((node->label())>0)){ // real tip node
			//std::ostringstream label_strm;
			//label_strm<<node->label();
			//std::ostringstream bl_strm;
			//bl_strm<< (node->parent_height() - node->height() + bl_above_parent)/4/npop;
			//return label_strm.str()+":"+bl_strm.str();
		//}
		//else if (node->is_root()){ // check if this is the root
				//return "("+ writeTree(node->first_child(),npop,0) +","+ writeTree(node->second_child(),npop,0)+")";
				//}
		//else{ // this is an interior node, but need to check if it is real, i.e. any of its children is a local
			//if (node->first_child()->local() && node->second_child()->local()){ // both children are local
				//std::ostringstream bl_strm;
				//bl_strm<< (node->parent_height() - node->height() + bl_above_parent)/4/npop;
				//return "("+ writeTree(node->first_child(),npop,0) +","+ writeTree(node->second_child(),npop,0)+"):"+bl_strm.str();;
			//}
			//else if(node->first_child()->local() && !node->second_child()->local()) { // first child is local, 
				//double local_bl_above_parent = bl_above_parent + node->parent_height() - node->height();
				//return writeTree(node->first_child(),npop, local_bl_above_parent);
			//}
			//else {// !node->first_child()->local() && node->second_child()->local() // second child is local, 
				//double local_bl_above_parent = bl_above_parent + node->parent_height() - node->height();
				//return writeTree(node->second_child(),npop, local_bl_above_parent);
			//}
			
		//}
	//}
//}

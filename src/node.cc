#include "node.h"
//#include<iostream>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<iostream>


Node::Node() { init(); };
Node::Node(double height) { init(height, true); };
Node::Node(double height, bool local) { init(height, local); };
Node::Node(double height, bool local, size_t last_update) { 
  init(height, local, last_update); };
Node::Node(double height, bool local, size_t last_update, size_t samples_below)
  { init(height, local, last_update, samples_below); };
Node::Node(double height, bool local, size_t last_update, size_t samples_below, double length_below)
  { init(height, local, last_update, samples_below, length_below); };
Node::Node(double height, bool local, size_t last_update, size_t samples_below, double length_below, int label)
  { init(height, local, last_update, samples_below, length_below, label); };
  
//Node::~Node() {delete[] marginal_likelihood;}
Node::~Node() {};

void Node::init(double height, bool local, size_t last_update,
                size_t samples_below, double length_below,int label) {
  this->set_height(height);
  this->set_local(local);
  this->set_parent(NULL);
  this->set_second_child(NULL);
  this->set_first_child(NULL);
  this->set_last_update(last_update);
  this->set_samples_below(samples_below);
  this->set_length_below(length_below);
  this->set_previous(NULL);
  this->set_next(NULL);
  this->set_label(label);
  this->mutation_state=false;
  this->marginal_likelihood[0]=0;
  this->marginal_likelihood[1]=0;
  this->set_mut_num(0);
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
  if ( this->height() == 0 ) return true;
  return false;
}






/**
 * Extract the string to represent the subtree that is descendant from node.

*/
std::string writeTree(Node * node, int npop){
	
	//if (node->tree_topo_bl.size()>0){
		//return node->tree_topo_bl;
	//}
	//else{
		if(node->first_child() == NULL && ((node->label())>0)){
		//if(node->label()>0){
			std::ostringstream label_strm;
			label_strm<<node->label();
			std::ostringstream bl_strm;
			bl_strm<< (node->parent_height() - node->height())/4/npop;
			node->tree_topo_bl = label_strm.str()+":"+bl_strm.str();
			return node->tree_topo_bl;
		}
		else{	
			//if (node->is_root()){
			if (!node->local()){
				node->tree_topo_bl="("+ writeTree(node->first_child(),npop) +","+ writeTree(node->second_child(),npop)+");";
				return node->tree_topo_bl;
				}
			else{
				std::ostringstream bl_strm;
				bl_strm<< (node->parent_height() - node->height())/4/npop;
				node->tree_topo_bl="("+ writeTree(node->first_child(),npop) +","+ writeTree(node->second_child(),npop)+"):"+bl_strm.str();;
				return node->tree_topo_bl;
			}
		}	
	//}
}

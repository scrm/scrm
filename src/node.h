/*
 * node.h
 *
 * A Node is the most elemental unit of a tree/forest. It represents the
 * point of coalescence of two branches, as well as the single branch above. 
 *
 */

#ifndef scrm_node
#define scrm_node

#include <cstddef>
#include <cfloat>
#include <stdexcept>
#include <iostream>
#include <cassert>
//#include <string>
#include <vector>
class Node
{
 public:                       

#ifdef UNITTEST
  friend class TestNode;
#endif
  size_t index; // this is the index of the node in the node container...
  
  Node();
  Node(double height);
  Node(double height, bool local);
  Node(double height, bool local, size_t last_update);
  Node(double height, bool local, size_t last_update, size_t samples_below);
  Node(double height, bool local, size_t last_update, size_t samples_below, double length_below);
  Node(double height, bool local, size_t last_update, size_t samples_below, double length_below, int label);

  ~Node();

  //Getters & Setters
  double height() const { return this->height_; }
  void set_height(const double &height) { this->height_ = height; }

  double parent_height() const {
    if ( this->is_root() ) return this->height();
    return this->parent()->height();
  }

  double height_above() const { return this->parent_height() - this->height(); }

  bool local() const { return this->local_; }
  void set_local(bool local) { this->local_ = local; }

  void make_local() { this->set_local(true); }
  void make_nonlocal(const size_t &current_base) { 
    if ( !local() ) return;
    set_last_update(current_base);
    this->set_local(false);
  }

	Node *parent() const; 
   void set_parent(Node *parent) { this->parent_ = parent; }
  //void set_parent_copy(Node *parent) { parent_ = parent; }


  Node *second_child() const { return this->second_child_; }
  void set_second_child(Node *second_child) { this->second_child_ = second_child; }

  Node *first_child() const { return this->first_child_; }
  void set_first_child(Node *first_child) { this->first_child_ = first_child; }

  size_t last_update() const { if ( local() ) return 0; 
    return(last_update_); }
  void set_last_update(size_t position) { this->last_update_ = position; }

  size_t samples_below() const { return samples_below_; }
  void set_samples_below(size_t samples) { samples_below_ = samples; }

  double length_below() const { return length_below_; }
  void set_length_below(double length) { length_below_ = length; }

  void change_child(Node* from, Node* to);
  int  numberOfChildren() const { 
    if (first_child() == NULL) return 0;
    else if (second_child() == NULL) return 1;
    else return 2;
  }
  
  void set_label(int label){label_=label;}
  int label() const{return label_;}

  bool is_root() const { return ( this->parent_ == NULL ); }
  bool in_sample() const;

  bool is_first() const { return( previous_ == NULL ); }
  bool is_last() const { return( next_ == NULL ); }

  void remove_child(Node* child);

  Node* next() const { 
    if ( next_ == NULL ) throw std::out_of_range("Node has no next node");
    return next_; 
  }
  Node* previous() const { 
    if ( previous_ == NULL ) throw std::out_of_range("Node has no previous node");
    return previous_; 
  }

  void set_next(Node* next) { next_ = next; }
  void set_previous(Node* previous) { previous_ = previous; }

  
  //std::string tree_topo_bl;
  bool mutation_state; // mutation state X = false denote homozygous site, X = true denote hetrozygous site
  double marginal_likelihood[2]; //marginal_likelihood[0] is the marginal probability of P(X = 0), 
								//marginal_likelihood[1] is the marginal probability of P(X = 1), 

  std::vector <int> descndnt;
  double mut_num() const { return this->mut_num_; }
  void set_mut_num(const double &mut_num) { this->mut_num_ = mut_num; }

 private:
 
  void init(double heigh=-1, 
            bool local=true, 
            size_t last_update = 0, 
            size_t samples_below=0, 
            double length_below=0,
            int label=0);
  int label_;
  double height_;        // The total height of the node
  bool   local_;        // Indicates if the branch above is local,
  // i.e. on the local tree
  size_t last_update_;   // The sequence position on which the branch above the node
  // was last checked for recombination events. Ignored for
  // local nodes, which are always up to date
  
  size_t samples_below_; // the number of sampled nodes in the subtree below this node
  double length_below_;  // the total length of local branches in the subtree below this node

  Node* next_;
  Node* previous_;

  //The tree structure
  Node *parent_;
  Node *first_child_;
  Node *second_child_;
  
  double mut_num_;
};
std::string writeTree(Node * node,int npop,double bl_above);

#endif

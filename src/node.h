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


class Node
{
 public:                       

#ifdef UNITTEST
  friend class TestNode;
#endif

  Node();
  Node(double height);
  Node(double height, bool local);
  Node(double height, bool local, size_t last_update);
  Node(double height, bool local, size_t last_update, size_t samples_below);
  Node(double height, bool local, size_t last_update, size_t samples_below, double length_below);
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
  void set_parent(Node *parent) { this->parent_ = parent; }; 

  Node *higher_child() const { return this->higher_child_; }
  void set_higher_child(Node *higher_child) { this->higher_child_ = higher_child; }

  Node *lower_child() const { return this->lower_child_; }
  void set_lower_child(Node *lower_child) { this->lower_child_ = lower_child; }

  size_t last_update() const { if ( local() ) return 0; 
    return(last_update_); }
  void set_last_update(size_t position) { this->last_update_ = position; }

  size_t samples_below() const { return samples_below_; }
  void set_samples_below(size_t samples) { samples_below_ = samples; }

  double length_below() const { return length_below_; }
  void set_length_below(double length) { length_below_ = length; }

  void change_child(Node* from, Node* to);
  void sort_children();
  int  numberOfChildren() const;

  bool is_root() const { return ( this->parent_ == NULL ); }
  bool in_sample() const;

  bool is_first() const { return( previous_ == NULL ); }
  bool is_last() const { return( next_ == NULL ); }

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

 private:
  void init(double heigh=-1, 
            bool local=true, 
            size_t last_update = 0, 
            size_t samples_below=0, 
            double length_below=0);

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
  Node *higher_child_;
  Node *lower_child_;   // If it has only one child, then it is the lower one
};

#endif

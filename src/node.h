/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013, 2014 Paul R. Staab, Sha (Joe) Zhu and Gerton Lunter
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

/*
 * node.h
 *
 * A Node is the most elemental unit of a tree/forest. It represents the
 * point of coalescence of two branches, as well as the single branch above. 
 *
 */

#ifndef scrm_node
#define scrm_node

#include "macros.h" // Needs to be before cassert

#include <cstddef>
#include <cfloat>
#include <stdexcept>
#include <iostream>
#include <cassert>

class Node
{
 public:                       

#ifdef UNITTEST
  friend class TestForest;
  friend class TestNode;
#endif
  
  Node();
  Node(double height);
  Node(double height, size_t label);

  ~Node();

  //Getters & Setters
  double height() const { return this->height_; }
  void set_height(const double height) { this->height_ = height; }

  double parent_height() const {
    if ( this->is_root() ) return this->height();
    return this->parent()->height();
  }

  double height_above() const { return this->parent_height() - this->height(); }

  size_t population() const { return population_; }
  void set_population(const size_t pop) { population_ = pop; }

  bool local() const { return (last_update_ == -1); }

  void make_local() { 
    last_update_ = -1; 
  }
  void make_nonlocal(const double current_base) { 
    assert( this->local() ); 
    set_last_update(current_base);
  }

	Node *parent() const {
    assert( this->parent_ != NULL ); 
    return this->parent_; 
  }
  void set_parent(Node *parent) { this->parent_ = parent; }

  Node *second_child() const { return this->second_child_; }
  void set_second_child(Node *second_child) { this->second_child_ = second_child; }

  Node *first_child() const { return this->first_child_; }
  void set_first_child(Node *first_child) { this->first_child_ = first_child; }

  double last_update() const { return last_update_; }

  double last_change() const { return last_change_; }
  void set_last_change(const double pos) { last_change_ = pos; }

  size_t samples_below() const { return samples_below_; }
  void set_samples_below(size_t samples) { samples_below_ = samples; }

  double length_below() const { return length_below_; }
  void set_length_below(const double length) { length_below_ = length; }

  void change_child(Node* from, Node* to);
  size_t countChildren(const bool only_local = false) const;
  
  void set_label(size_t label) { label_ = label; }
  size_t label() const { return label_; }

  void set_OF_label(size_t OF_label) { OF_label_ = OF_label; }
  size_t OF_label() const { return OF_label_; }
  
  bool is_root() const { return ( this->parent_ == NULL ); }
  bool in_sample() const {
    return ( this->label() != 0 ); 
  }

  bool is_migrating() const; 

  bool is_first() const { return( previous_ == NULL ); }
  bool is_last() const { return( next_ == NULL ); }

  // Uminportant Nodes are nodes that sit at the top of the single 
  // top branch of a tree after it got cut away from the primary tree. 
  // These Nodes can be reused or removed if they are involved in an event.
  bool is_unimportant() const {
    return (this->countChildren() == 1 && !this->is_migrating()); 
  }

  bool is_contemporary(const double time) {
    return ( time <= height() && height() <= parent_height() ); 
  }

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

  // Navigate on local tree
  Node *getLocalParent() const;
  Node *getLocalChild1() const;
  Node *getLocalChild2() const;

 private:
  void init(double heigh=-1, size_t label=0);
  void set_last_update(const double position) { last_update_ = position; }; 

  size_t label_;
  size_t OF_label_;
  double height_;        // The total height of the node
  double last_update_;   // The sequence position on which the branch above the node
                         // was last checked for recombination events or 0 if
                         // the node is local 

  double last_change_;

  size_t population_;    // The number of the population the node belong to. 
  
  size_t samples_below_; // the number of sampled nodes in the subtree below this node
  double length_below_;  // the total length of local branches in the subtree below this node

  Node* next_;
  Node* previous_;

  //The tree structure
  Node *parent_;
  Node *first_child_;
  Node *second_child_;
};

inline bool Node::is_migrating() const { 
  if ( this->countChildren() != 1 ) return false;
  return ( this->population() != this->first_child()->population() );
} 

inline size_t Node::countChildren(const bool only_local) const { 
  if (first_child() == NULL) return 0;
  if (!only_local) {
    if (second_child() == NULL) return 1;
    else return 2;
  } else {
    if (second_child() == NULL) {
      if (first_child()->local()) return 1; 
      else return 0;
    } else {
      return first_child()->local() + second_child()->local();
    }
  }
}

#endif

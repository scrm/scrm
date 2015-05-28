/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013, 2014 Paul R. Staab, Sha (Joe) Zhu, Dirk Metzler and Gerton Lunter
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
#include<sstream>
#include<fstream>
#include<iomanip>
#include<iostream>


Node::Node() { init(); }
Node::Node(double height) { init(height); }
Node::Node(double height, size_t label) { init(height, label); }
Node::~Node() {}

void Node::init(double height, size_t label) {
  this->set_last_update(0);
  this->set_population(0);

  this->set_height(height);
  this->set_label(label);
  if (label == 0) this->set_samples_below(1);
  else this->set_samples_below(0); 
  this->set_length_below(0);
  this->set_last_change(0);

  this->set_parent(NULL);
  this->set_second_child(NULL);
  this->set_first_child(NULL);
  this->set_previous(NULL);
  this->set_next(NULL);
}
  

void Node::change_child(Node* from, Node* to) {
  if ( this->first_child() == from ) {
    if (to != NULL) this->set_first_child(to);
    else {
      set_first_child(second_child());
      set_second_child(NULL);
    }
  }
  else if ( this->second_child() == from )
    this->set_second_child(to);
  else {
    dout << "Error when changing child of " << this << " form "
         << from << " to " << to << std::endl;
    dout << "Children are " << this->first_child() << " and "
         << this->second_child() << std::endl;
    throw std::invalid_argument("Can't find child node to replace");
  }
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

/**
 * @brief Returns is the parent of this node on the local tree.
 *
 * This should only be executed on local nodes!
 *
 * @return The node that is this node's parent on the local tree. 
 */
Node* Node::getLocalParent() const {
  assert( this->local() );
  assert( !this->is_root());
  assert( this->parent()->countChildren() > 0 );
  if (parent()->countChildren(true) == 2) return parent(); 
  return parent()->getLocalParent();
}

/**
 * @brief Returns one child of this node when looking 
 * only at the local tree.
 * 
 * The up to two children of a node are in basically 
 * arbitrary order. The only difference between child_1 and
 * child_2 is that if a node has only one child, than
 * it is always child_1.
 *
 * This should only be executed on local nodes!
 *
 * @return NULL if this node has no children on the local tree,
 * or a pointer to the child otherwise. 
 */
Node* Node::getLocalChild1() const {
  if (first_child() == NULL || !first_child()->local()) return NULL;
  if (first_child()->countChildren(true) == 1) {
    if (first_child()->first_child()->local()) return first_child()->getLocalChild1();
    else return first_child()->getLocalChild2();
  }
  // Child has 0 or 2 local children => it is part of local tree 
  return first_child();
}

/**
 * @brief Returns one child of this node when looking 
 * only at the local tree.
 * 
 * The up to two children of a node are in basically 
 * arbitrary order. The only difference between child_1 and
 * child_2 is that if a node has only one child, than
 * it is always child_1.
 *
 * This should only be executed on local nodes!
 *
 * @return NULL if this node has less than two children on the local tree,
 * or a pointer to the child otherwise. 
 */
Node* Node::getLocalChild2() const {
  if (second_child() == NULL || !second_child()->local()) return NULL;
  if (second_child()->countChildren(true) == 1) {
    if (second_child()->first_child()->local()) return second_child()->getLocalChild1();
    else return second_child()->getLocalChild2();
  }
  // Child has 0 or 2 local children => it is part of local tree 
  return second_child();
}


void Node::extract_bl_and_label ( std::string::iterator in_it ){
  // Going backwards, extract branch length first, then the node label
  std::string::iterator bl_start = in_it;
  //for (; (*(bl_start-1)) != ':'; --bl_start ){ }
  while ((*(bl_start-1)) != ':')
    --bl_start;
  double bl = strtod( &(*bl_start), NULL );
  this->set_bl ( bl );

  std::string::iterator label_start = (bl_start-2);

  while ( (*(label_start)) != ',' &&  ( *(label_start)) != '(' &&  (*(label_start)) != ')' )
    --label_start;

  this->set_label ( ( (*(label_start)) == ')' ? 0 /*! Label internal nodes */
                                                : strtol ( &(*(label_start+1)), NULL , 10)) ); /*! Label tip nodes */
}


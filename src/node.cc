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
  this->set_last_update(-1);
  this->set_population(0);

  this->set_height(height);
  this->set_label(label);
  if (label == 0) this->set_samples_below(1);
  else this->set_samples_below(0); 
  this->set_OF_label(label);
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

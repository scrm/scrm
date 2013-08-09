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

#ifndef scrm_src_event
#define scrm_src_event

#include "node.h"

class Event {
 public:
  Event() {
    type_  = 0;
    time_  = -1;
    node_ = NULL;
    mig_pop_ = -1;
  };

  Event(double time) {
    type_ = 0;
    time_ = time;
    node_ = NULL;
    mig_pop_ = -1;
  }
  
  double time() { return time_; }
  Node* node() { return node_; }
  size_t mig_pop() { return mig_pop_; }
  size_t type() { return type_; }

  bool isNoEvent() { return (type_ == 0); }
  bool isCoalescence() { return (type_ == 1); }
  bool isPwCoalescence() { return (type_ == 2); }
  bool isMigration() { return (type_ == 3); }
  bool isRecombination() { return (type_ == 4); }
  
  void set_time(const double &time) { time_ = time; }

  Event setToCoalescence(Node* node) {
    type_  = 1;
    node_ = node;
    return *this;
  }
  Event setToPwCoalescence() {
    type_  = 2;
    return *this;
  }
  Event setToMigration(Node* node, const size_t &mig_pop) {
    type_  = 3;
    node_ = node;
    mig_pop_ = mig_pop;
    return *this;
  }
  Event setToRecombination(Node* node) {
    type_  = 4;
    node_ = node;
    return *this;
  }

 private:
  size_t type_;
  double time_;
  Node* node_;
  size_t mig_pop_;
};


#endif

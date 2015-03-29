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

#ifndef scrm_src_event
#define scrm_src_event

#include <iostream>
#include "node.h"

class Event {
 public:
  Event() {
    type_  = 0;
    time_  = -1;
    node_  = NULL;
    active_node_nr_ = -1;
    mig_pop_ = -1;
  };

  Event(double time) {
    type_ = 0;
    time_ = time;
    node_ = NULL;
    mig_pop_ = -1;
    active_node_nr_ = -1;
  }
  
  friend std::ostream& operator<< (std::ostream& stream, const Event& event);

  double time() const { return time_; }
  Node* node() const { return node_; }
  size_t mig_pop() const { return mig_pop_; }
  size_t type() const { return type_; }
  size_t active_node_nr() const { return active_node_nr_; }

  bool isNoEvent() const { return (type_ == 0); }
  bool isCoalescence() const { return (type_ == 1); }
  bool isPwCoalescence() const { return (type_ == 2); }
  bool isMigration() const { return (type_ == 3); }
  bool isRecombination() const { return (type_ == 4); }
  
  void set_time(const double time) { time_ = time; }

  void setToCoalescence(Node *node, const size_t active_node_nr) {
    type_  = 1;
    node_  = node;
    active_node_nr_ = active_node_nr;
  }
  void setToPwCoalescence() {
    type_  = 2;
  }
  void setToMigration(Node *node, const size_t active_node_nr, const size_t mig_pop) {
    type_  = 3;
    node_  = node;
    active_node_nr_ = active_node_nr;
    mig_pop_ = mig_pop;
  }
  void setToRecombination(Node *node, const size_t active_node_nr) {
    type_  = 4;
    node_  = node;
    active_node_nr_ = active_node_nr;
  }

 private:
  size_t type_;
  size_t active_node_nr_;
  double time_;
  size_t mig_pop_;
  Node* node_;
};

inline std::ostream& operator<< (std::ostream& stream, const Event& event) {
  if (event.isNoEvent()) {
    stream << "No Event";
    return stream;
  }

  stream << "Event at time " << event.time() << ": ";
  if (event.isCoalescence()) stream << "Coalesence of a" << event.active_node_nr();
  else if (event.isPwCoalescence()) stream << "Pair-wise coalescence";
  else if (event.isMigration()) stream << "Migration of Node " << event.node() 
                                        << " into pop " << event.mig_pop();
  else if (event.isRecombination()) stream << "Recombination of Node " << event.node();
  return stream;
}

#endif

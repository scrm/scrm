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

#ifndef scrm_src_node_container
#define scrm_src_node_container

#include "macros.h" // Needs to be before cassert

#include <vector>
#include <stack>
#include <stdexcept>
#include <cfloat>
#include <cassert>
#include <iostream>
#include <map>
#include <algorithm>
#include <list>

#include "node.h"

class NodeIterator;
class ConstNodeIterator;
class ReverseConstNodeIterator;

class NodeContainer {
 public:
  NodeContainer();
  ~NodeContainer() {
    clear();
    for (std::vector<Node>* lane : node_lanes_) delete lane;
  };

  NodeContainer& operator=(NodeContainer nc) {
    swap(*this, nc);
    return(*this);
  };

  NodeContainer(const NodeContainer &nc);

  NodeIterator iterator();
  NodeIterator iterator(Node* node);
  ConstNodeIterator iterator() const;
  ConstNodeIterator iterator(Node* node) const;
  ReverseConstNodeIterator reverse_iterator() const;
  ReverseConstNodeIterator reverse_iterator(Node* node) const;

  // Create Nodes
  Node* createNode(double height, size_t label = 0) {
    // Use the slot of a previously deleted node if possible
    if (free_slots_.size() > 0) {
      Node* node = free_slots_.top();
      free_slots_.pop();
      *node = Node(height, label);
      return node;
    }

    // Otherwise, use a new slot
    if (node_counter_ >= 10000) {
      ++lane_counter_;
      node_counter_ = 0;
      if (lane_counter_ == node_lanes_.size()) {
        std::vector<Node>* new_lane = new std::vector<Node>();
        new_lane->reserve(10000);
        node_lanes_.push_back(new_lane);
      }
    }
    ++node_counter_;
    node_lanes_.at(lane_counter_)->push_back(Node(height, label));
    return &*(node_lanes_.at(lane_counter_)->end() - 1);
  }

  // Create Nodes
  Node* createNode(const Node copiedNode) {
    // Use the slot of a previously deleted node if possible
    if (free_slots_.size() > 0) {
      Node* node = free_slots_.top();
      free_slots_.pop();
      *node = Node(copiedNode);
      return node;
    }

    // Otherwise, use a new slot
    if (node_counter_ >= 10000) {
      ++lane_counter_;
      node_counter_ = 0;
      if (lane_counter_ == node_lanes_.size()) {
        std::vector<Node>* new_lane = new std::vector<Node>();
        new_lane->reserve(10000);
        node_lanes_.push_back(new_lane);
      }
    }
    ++node_counter_;
    node_lanes_.at(lane_counter_)->push_back(copiedNode);
    return &*(node_lanes_.at(lane_counter_)->end() - 1);
  }

  void push_back(Node* node);
  void push_front(Node* node);
  void add(Node* node, Node* after_node=NULL);
  void remove(Node *node, const bool &del=true);
  void move(Node *node, const double new_height);
  void clear();

  Node* at(size_t nr) const;
  Node const* get(size_t nr) const { return at(nr); };

  Node* first() const { return first_node_; };
  Node* last() const { return last_node_; };

  size_t size() const { return size_; };
  bool sorted() const;

#ifdef UNITTEST
  friend class TestNodeContainer;
#endif
  friend class NodeIterator;
  friend class ConstNodeIterator;
  friend class ReverseConstNodeIterator;
  friend std::ostream& operator<< (std::ostream& stream, const NodeContainer& nc);

 private:
  friend void swap(NodeContainer& first, NodeContainer& second);
  //const std::vector<Node*> nodes() const { return nodes_; }

  void add_before(Node* add, Node* next_node);

  Node* first_node_;
  Node* last_node_;

  void set_first(Node* node) { first_node_ = node; }
  void set_last(Node* node) { last_node_ = node; }

  Node* unsorted_node_;
  size_t size_;

  // Storing the nodes in lanes a 10k nodes
  std::vector<std::vector<Node>*> node_lanes_;
  std::stack<Node*> free_slots_;
  size_t node_counter_;
  size_t lane_counter_;
};


class NodeIterator {
 public:
  NodeIterator() { current_node_ = NULL; };
  NodeIterator(NodeContainer& nc) { current_node_ = nc.first(); };
  NodeIterator(Node* node) { current_node_ = node; };
  ~NodeIterator() {};

  Node* operator*() {
    if (current_node_ == NULL)
      throw std::out_of_range("Node iterator out of range");
    return current_node_;
  }

  Node* operator++() {
    if (current_node_ == NULL)
      throw std::out_of_range("Node iterator out of range");

    if ( current_node_->is_last() ) current_node_ = NULL;
    else current_node_ = current_node_->next();
    return current_node_;
  }

  Node* operator--() {
    if (current_node_ == NULL)
      throw std::out_of_range("Node iterator out of range");

    if ( current_node_->is_first() ) current_node_ = NULL;
    else current_node_ = current_node_->previous();
    return current_node_;
  }

  Node* operator++(int) {
    if (current_node_ == NULL)
      throw std::out_of_range("Node iterator out of range");

    Node* ret = current_node_;
    if ( current_node_->is_last() ) current_node_ = NULL;
    else current_node_ = current_node_->next();
    return ret;
  }

  bool good() const {
    return current_node_ != NULL;
  }

  double height() const {
    if ( good() ) return current_node_->height();
    else return DBL_MAX;
  }

#ifdef UNITTEST
  friend class TestNodeContainer;
#endif

 private:
  Node* current_node_;
};


class ConstNodeIterator {
 public:
  ConstNodeIterator() { current_node_ = NULL; };
  ConstNodeIterator(const NodeContainer& nc) { current_node_ = nc.first(); };
  ConstNodeIterator(Node const* node) { current_node_ = node; };
  ~ConstNodeIterator() {};

  Node const* operator*() {
    if (current_node_ == NULL)
      throw std::out_of_range("Node iterator out of range");
    return current_node_;
  }

  Node const* operator++() {
    if (current_node_ == NULL)
      throw std::out_of_range("Node iterator out of range");

    if ( current_node_->is_last() ) current_node_ = NULL;
    else current_node_ = current_node_->next();

    return current_node_;
  }

  Node const* operator--() {
    if (current_node_ == NULL)
      throw std::out_of_range("Node iterator out of range");

    if ( current_node_->is_first() ) current_node_ = NULL;
    else current_node_ = current_node_->previous();
    return current_node_;
  }

  Node const* operator++(int) {
    if (current_node_ == NULL)
      throw std::out_of_range("Node iterator out of range");

    Node const* ret = current_node_;
    if ( current_node_->is_last() ) current_node_ = NULL;
    else current_node_ = current_node_->next();
    return ret;
  }

  bool good() const {
    return current_node_ != NULL;
  }

  double height() const {
    if ( good() ) return current_node_->height();
    else return DBL_MAX;
  }

 private:
  Node const* current_node_;
};


class ReverseConstNodeIterator {
 public:
  ReverseConstNodeIterator() { current_node_ = NULL; };
  ReverseConstNodeIterator(const NodeContainer &nc) {  current_node_ = nc.last(); };
  ReverseConstNodeIterator(Node const* node) { current_node_ = node; };
  ~ReverseConstNodeIterator() {};

  Node const* operator*() {
    if (current_node_ == NULL)
      throw std::out_of_range("Node iterator out of range");
    return current_node_;
  }

  Node const* operator++() {
    if (current_node_ == NULL)
      throw std::out_of_range("Node iterator out of range");

    if ( current_node_->is_first() ) current_node_ = NULL;
    else current_node_ = current_node_->previous();
    return current_node_;
  }

  Node const* operator++(int) {
    if (current_node_ == NULL)
      throw std::out_of_range("Node iterator out of range");

    Node const* ret = current_node_;
    if ( current_node_->is_first() ) current_node_ = NULL;
    else current_node_ = current_node_->previous();
    return ret;
  }

  bool good() const {
    return current_node_ != NULL;
  }

 private:
  Node const* current_node_;
};


inline NodeIterator NodeContainer::iterator() { return NodeIterator(*this); }
inline NodeIterator NodeContainer::iterator(Node* node) { return NodeIterator(node); }
inline ConstNodeIterator NodeContainer::iterator() const { return ConstNodeIterator(*this); }
inline ConstNodeIterator NodeContainer::iterator(Node* node) const { return ConstNodeIterator(node); }
inline ReverseConstNodeIterator NodeContainer::reverse_iterator() const { return ReverseConstNodeIterator(*this); }
inline ReverseConstNodeIterator NodeContainer::reverse_iterator(Node* node) const {return ReverseConstNodeIterator(node); }
#endif

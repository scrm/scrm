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

#include "node_container.h"

/*******************************************************
 * Con- & Destructor
 *******************************************************/

NodeContainer::NodeContainer() {
  set_first(NULL);
  set_last(NULL);
  unsorted_node_ = NULL;
  size_ = 0;

  node_counter_ = 0;
  lane_counter_ = 0; 
  std::vector<Node>* new_lane = new std::vector<Node>();
  new_lane->reserve(10000);
  node_lanes_.push_back(new_lane);
}

NodeContainer::NodeContainer(const NodeContainer &nc) {
  size_ = 0;
  set_first(NULL);
  set_last(NULL);
  this->unsorted_node_ = NULL;
  this->last_node_ = NULL;
  this->first_node_ = NULL;

  std::map<Node const*, Node*> node_mapping;  
  node_mapping[NULL] = NULL;

  for (auto it = nc.iterator(); it.good(); ++it) {
    Node *node = new Node(**it);
    node_mapping[*it] = node;
    add(node);
  }
  assert( this->sorted() );

  for (auto it = iterator(); it.good(); ++it) {
    if (!(*it)->is_root()) (*it)->set_parent(node_mapping[(*it)->parent()]);
    (*it)->set_first_child(node_mapping[(*it)->first_child()]);
    (*it)->set_second_child(node_mapping[(*it)->second_child()]);
  }

  unsorted_node_ = node_mapping[nc.unsorted_node_];
}

/*******************************************************
 * Management of Nodes
 *******************************************************/

Node* NodeContainer::at(size_t nr) const {
  Node* current = first();

  for (size_t i=0; i < nr; ++i) {
    assert(current != NULL);
    current = current->next();
  }

  if ( current == NULL ) throw std::out_of_range("NodeContainer out of range"); 
  return current;
}

// Adds 'node' to the container
void NodeContainer::push_back( Node* node ) {
    ++size_;
    if ( this->first() == NULL ) {
        this->set_first( node );
        this->set_last( node );
        return;
    }
    assert( this->first() != NULL );

    //Adding to the End, similar to vector::push_back
    node->set_previous(this->last());
    node->set_next(NULL);
    this->last()->set_next(node);
    this->set_last(node);
    return;
}

// Adds 'node' to the container
void NodeContainer::push_front( Node* node ) {
    ++size_;
    if ( this->first() == NULL ) {
        this->set_first( node );
        this->set_last( node );
        return;
    }
    assert( this->first() != NULL );

    //Adding to the End, similar to vector::push_back
    node->set_next(first());
    node->set_previous(NULL);
    first()->set_previous(node);
    this->set_first(node);
    return;
}



// Adds 'node' to the container
// If you know that the node is higher in the tree than a other node,
// than you can specify the latter as 'after_node' to speedup the process.
void NodeContainer::add(Node* node, Node* after_node) {
  ++size_;

  if (first() == NULL) {
    this->set_first(node);
    this->set_last(node);
    return;
  }
  assert(first() != NULL);

  // Before first node?
  if (node->height() <= first()->height()) {
    node->set_next(first());
    node->set_previous(NULL);
    first()->set_previous(node);
    this->set_first(node);
    assert( this->sorted() );
    return;
  }

  // After last node?
  if ( node->height() >= last()->height() ) {
    node->set_previous(last());
    node->set_next(NULL);
    last()->set_next(node);
    this->set_last(node);
    assert( this->sorted() );
    return;
  }

  assert( after_node == NULL || node->height() >= after_node->height() );


  if (after_node == NULL) after_node = first();
  Node* current = after_node;
  // Find position in between
  while ( current->height() <= node->height() ) {
    assert( !current->is_last() );
    if ( !current->is_root() ) {
      if ( current->parent_height() < node->height() ) { 
        current = current->parent();
        continue;
      }
    }
    current = current->next();
  }
 
  // And add the node;
  this->add_before(node, current);
  assert( this->sorted() );
}


void NodeContainer::remove(Node *node, const bool &del) {
  --size_; 
  if ( node->is_first() && node->is_last() ) {
    this->set_first(NULL);
    this->set_last(NULL);
  }
  else if ( node->is_first() ) {
    this->set_first(node->next());
    node->next()->set_previous(NULL);
  }
  else if ( node->is_last() ) {
    this->set_last(node->previous());
    node->previous()->set_next(NULL);
  }
  else {
    node->previous()->set_next(node->next());
    node->next()->set_previous(node->previous());
  }
  
  if (del) free_slots_.push(node);
  assert( this->sorted() );
}


void NodeContainer::move(Node *node, const double new_height) {
  assert( node != NULL );

  // Stupid edge case first: We may have only one node.
  if ( node->is_first() && node->is_last() ) {
    node->set_height(new_height);
    return;
  }

  // Remove from old place
  remove(node, false);

  // Add at new place
  Node* current = NULL;
  if ( node->height() < new_height ) {
    if ( node->is_first() ) current = NULL;
    else current = node->previous();
  } else {
    current = first();
  }

  node->set_height(new_height);
  this->add(node, current);
  assert( this->sorted() );
}


void NodeContainer::clear() {
  set_first(NULL);
  set_last(NULL);
  this->size_ = 0;
  this->node_counter_ = 0;
  this->lane_counter_ = 0;

  // Clear free_slots_
  std::stack<Node*>().swap(free_slots_);
}

void NodeContainer::add_before(Node* add, Node* next_node){
  add->set_next(next_node);
  add->set_previous(next_node->previous());

  if ( add->previous() != NULL ) add->previous()->set_next(add);
  next_node->set_previous(add);
  if ( add->is_last() ) this->set_last(add);
}



/*******************************************************
 * Consistency checking
 *******************************************************/

bool NodeContainer::sorted() const {
  Node* current = first();
  if ( !current->is_first() ) {
    dout << "NodeContainer: First Node is not first" << std::endl;
    return 0;
  }

  while ( !current->is_last() ) {
    current = current->next();
    if ( current->height() < current->previous()->height() ) {
      dout << "NodeContainer: Nodes not sorted" << std::endl;
      return 0;
    }
    if ( current == current->previous() ) {
      dout << "NodeContainer: Fatal loop detected" << std::endl;
      return 0;
    }
  }
  
  if ( !current->is_last() ) {
    dout << "NodeContainer: Last Node not last" << std::endl;
    return 0;
  }

  return 1;
}


void swap(NodeContainer& first, NodeContainer& second) {
  using std::swap;
  swap(first.first_node_, second.first_node_);
  swap(first.last_node_, second.last_node_);
  swap(first.size_, second.size_);
  swap(first.unsorted_node_, second.unsorted_node_);
}


std::ostream& operator<< (std::ostream& stream, const NodeContainer& nc) {
  for ( ConstNodeIterator it = nc.iterator(); it.good(); ++it ) {
    stream << *it << "(" << (*it)->height() << ")";
    if (*it != nc.last()) stream << " <--> "; 
  }
  return stream;
}


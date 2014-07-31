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

#ifndef scrm_src_contemporaries_container
#define scrm_src_contemporaries_container

#include <vector>
#include <unordered_set>
#include <cassert>

#include "node.h"
#include "random/random_generator.h"

#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

class ContemporariesIterator {
 public:
  ContemporariesIterator(std::unordered_set<Node*>::iterator it) {
    it_set_ = it; 
  };

  Node* operator*() const { return *it_set_; }
  Node* operator++() { ++it_set_; return *it_set_; }
  Node* operator++(int) { Node* node = *it_set_; ++it_set_; return node; }
  
  bool operator==(const ContemporariesIterator &other) const {
    return (it_set_ == other.it_set_); 
  }

  bool operator!=(const ContemporariesIterator &other) const {
    return !(*this == other);
  }

 private:
  std::unordered_set<Node*>::iterator it_set_;
};

class ContemporariesConstIterator {
 public:
  ContemporariesConstIterator(std::unordered_set<Node*>::const_iterator it) {
    it_set_ = it; 
  };

  Node const* operator*() const { return *it_set_; }
  Node const* operator++() { ++it_set_; return *it_set_; }
  Node const* operator++(int) { Node* node = *it_set_; ++it_set_; return node; }
  
  bool operator==(const ContemporariesConstIterator &other) const {
    return (it_set_ == other.it_set_); 
  }

  bool operator!=(const ContemporariesConstIterator &other) const {
    return !(*this == other);
  }

 private:
  std::unordered_set<Node*>::const_iterator it_set_;
};

class ContemporariesContainer {
 public:
  ContemporariesContainer();
  ContemporariesContainer(const size_t pop_number, 
                          const size_t sample_number,
                          RandomGenerator *rg);

  void add(Node* node);
  void remove(Node *node);
  void replaceChildren(Node *add_node);
  void replace(Node *add_node, Node *del_node_1, Node *del_node_2 = NULL);
  void clear();
  void buffer(const double current_time);


  Node* sample(const size_t pop) const;
  size_t size(const size_t pop) const;
  double buffer_time() const { return buffer_time_; }
  bool empty() const;

  ContemporariesConstIterator begin(const size_t pop) const { 
    return ContemporariesConstIterator(contemporaries_set().at(pop).cbegin()); 
  }
  ContemporariesConstIterator end(const size_t pop) const { 
    return ContemporariesConstIterator(contemporaries_set().at(pop).cend()); 
  }

  ContemporariesConstIterator buffer_begin(const size_t pop) const { 
    return ContemporariesConstIterator(buffer_set().at(pop).cbegin()); 
  }
  ContemporariesConstIterator buffer_end(const size_t pop) const { 
    return ContemporariesConstIterator(buffer_set().at(pop).cend()); 
  }
  ContemporariesIterator buffer_begin(const size_t pop) { 
    return ContemporariesIterator(buffer_set().at(pop).begin()); 
  }
  ContemporariesIterator buffer_end(const size_t pop) { 
    return ContemporariesIterator(buffer_set().at(pop).end()); 
  }

 private:
  std::vector<std::unordered_set<Node*> > &contemporaries_set() {
    if (use_first_) return contemporaries_set1_;
    else return contemporaries_set2_;
  }
  const std::vector<std::unordered_set<Node*> > &contemporaries_set() const {
    if (use_first_) return contemporaries_set1_;
    else return contemporaries_set2_;
  }
  std::vector<std::unordered_set<Node*> > &buffer_set() {
    if (!use_first_) return contemporaries_set1_;
    else return contemporaries_set2_;
  }
  const std::vector<std::unordered_set<Node*> > &buffer_set() const {
    if (!use_first_) return contemporaries_set1_;
    else return contemporaries_set2_;
  }

  std::vector<std::unordered_set<Node*> > contemporaries_set1_;
  std::vector<std::unordered_set<Node*> > contemporaries_set2_;
  bool use_first_;
  double buffer_time_;
  RandomGenerator* rg_;
};

inline ContemporariesContainer::ContemporariesContainer() {
  contemporaries_set1_ = std::vector<std::unordered_set<Node*> >(1, std::unordered_set<Node*>(100));
  contemporaries_set2_ = std::vector<std::unordered_set<Node*> >(1, std::unordered_set<Node*>(100));
  rg_ = NULL;
  use_first_ = true;
  buffer_time_ = -1;
}

inline ContemporariesContainer::ContemporariesContainer(const size_t pop_number, 
                                                        const size_t sample_number,
                                                        RandomGenerator* rg) {
  size_t bucket_nr = std::ceil((sample_number + 200) * 1.4);
  contemporaries_set1_ = std::vector<std::unordered_set<Node*> >(pop_number, std::unordered_set<Node*>(bucket_nr));
  contemporaries_set2_ = std::vector<std::unordered_set<Node*> >(pop_number, std::unordered_set<Node*>(bucket_nr));
  this->rg_ = rg;
  use_first_ = true;
  buffer_time_ = DBL_MAX;
}

inline void ContemporariesContainer::add(Node* node) {
  assert(node != NULL);
  assert(!node->is_root());
  contemporaries_set().at(node->population()).insert(node);
};

inline void ContemporariesContainer::remove(Node* node) {
  assert(node != NULL);
  contemporaries_set().at(node->population()).erase(node);
}

inline void ContemporariesContainer::replaceChildren(Node *add_node) {
  replace(add_node, add_node->first_child(), add_node->second_child());  
};

inline void ContemporariesContainer::replace(Node *add_node, Node *del_node_1, Node *del_node_2) {
  assert(add_node != NULL);
  if (del_node_1 != NULL) remove(del_node_1);
  if (del_node_2 != NULL) remove(del_node_2);
  if (!add_node->is_root()) add(add_node);
};

inline void ContemporariesContainer::clear() {
  for (auto it = contemporaries_set().begin(); it != contemporaries_set().end(); ++it) {
    if ((*it).size() > 0) (*it).clear();
  }
}

inline size_t ContemporariesContainer::size(const size_t pop) const {
  return contemporaries_set().at(pop).size();
};

/**
 * @brief Function that buffers the current state of contemporaries for later
 * use. 
 *
 * Be careful not to change the tree at the time of buffer, as this
 * will not be reflected in the buffered contemporaries.
 *
 * @param current_time The time for with the current state is valid.
 *
 * @return 
 */
inline void ContemporariesContainer::buffer(const double current_time) {
  buffer_time_ = current_time;
  use_first_ = 1 - use_first_; 
  this->clear();
};

// Uniformly samples a random node from the current contemporaries.
// Distribution checked.
inline Node* ContemporariesContainer::sample(const size_t pop) const {
  assert( pop < contemporaries_set().size() );
  assert( this->size(pop) > 0 );

  // Sample the position of the Node we return
  size_t sample = rg_->sampleInt(this->size(pop));
  
  for (auto it = contemporaries_set().at(pop).begin(); 
       it != contemporaries_set().at(pop).end(); ++it) {
    assert( *it != NULL );
    if ( sample == 0 ) return (*it);
    --sample;
  }

  throw std::logic_error("Failed to find the contemporary I wanted to sample.");
  return NULL;
}

inline bool ContemporariesContainer::empty() const {
  for (auto pop_container : contemporaries_set()) {
    if ( !pop_container.empty() ) return false;
  }
  return true;
};
#endif

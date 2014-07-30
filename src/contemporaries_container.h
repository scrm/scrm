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
  ContemporariesIterator(std::unordered_set<Node*>::const_iterator it) {
    it_set_ = it; 
  };

  Node const* operator*() const { return *it_set_; }
  Node const* operator++() { ++it_set_; return *it_set_; }
  Node const* operator++(int) { Node* node = *it_set_; ++it_set_; return node; }
  
  bool operator==(const ContemporariesIterator &other) const {
    return (it_set_ == other.it_set_); 
  }

  bool operator!=(const ContemporariesIterator &other) const {
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
  void replace(Node *del_node, Node *add_node);
  void clear();

  Node* sample(const size_t pop) const;
  size_t size(const size_t pop) const;

  ContemporariesIterator begin(const size_t pop) const { 
    return ContemporariesIterator(contemporaries_set_.at(pop).cbegin()); 
  }
  ContemporariesIterator end(const size_t pop) const { 
    return ContemporariesIterator(contemporaries_set_.at(pop).cend()); 
  }

 private:
  std::vector<std::unordered_set<Node*> > contemporaries_set_;
  RandomGenerator* rg_;
};

inline ContemporariesContainer::ContemporariesContainer() {
  contemporaries_set_ = std::vector<std::unordered_set<Node*> >(1, std::unordered_set<Node*>(100));
  rg_ = NULL;
}

inline ContemporariesContainer::ContemporariesContainer(const size_t pop_number, 
                                                        const size_t sample_number,
                                                        RandomGenerator* rg) {
  size_t bucket_nr = std::ceil((sample_number + 200) * 1.4);
  contemporaries_set_ = std::vector<std::unordered_set<Node*> >(pop_number, std::unordered_set<Node*>(bucket_nr));
  this->rg_ = rg;
}


inline void ContemporariesContainer::add(Node* node) {
  assert(node != NULL);
  contemporaries_set_.at(node->population()).insert(node);
};

inline void ContemporariesContainer::remove(Node* node) {
  assert(node != NULL);
  contemporaries_set_.at(node->population()).erase(node);
}

inline void ContemporariesContainer::clear() {
  for (auto it = contemporaries_set_.begin(); it != contemporaries_set_.end(); ++it) {
    if ((*it).size() > 0) (*it).clear();
  }
}

inline size_t ContemporariesContainer::size(const size_t pop) const {
  return contemporaries_set_.at(pop).size();
};

// Uniformly samples a random node from the current contemporaries.
// Distribution checked.
inline Node* ContemporariesContainer::sample(const size_t pop) const {
  assert( pop < contemporaries_set_.size() );
  assert( this->size(pop) > 0 );

  // Sample the position of the Node we return
  size_t sample = rg_->sampleInt(this->size(pop));
  
  for (auto it = contemporaries_set_.at(pop).begin(); 
       it != contemporaries_set_.at(pop).end(); ++it) {
    assert( *it != NULL );
    if ( sample == 0 ) return (*it);
    --sample;
  }

  throw std::logic_error("Failed to find the contemporary I wanted to sample.");
  return NULL;
}

#endif

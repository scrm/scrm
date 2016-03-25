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

#ifndef scrm_src_contemporaries_container
#define scrm_src_contemporaries_container

#include "macros.h" // Needs to be before cassert

#include <vector>
#include <unordered_set>
#include <cassert>
#include <algorithm>

#include "node.h"
#include "random/random_generator.h"

class ContemporariesIterator {
 public:
  ContemporariesIterator(std::unordered_set<Node*>::iterator it) : it_set_(it), use_set_(true) { };
  ContemporariesIterator(std::vector<Node*>::iterator it) : it_vec_(it), use_set_(false) { };

  Node* operator*() const {
    if (use_set_) return *it_set_;
    else return *it_vec_;
  }

  ContemporariesIterator& operator++() {
    if (use_set_) (++it_set_);
    else ++it_vec_;
    return *this;
  }

  bool operator==(const ContemporariesIterator &other) const {
    if (use_set_) return (it_set_ == other.it_set_);
    else return (it_vec_ == other.it_vec_);
  }

  bool operator!=(const ContemporariesIterator &other) const {
    return !(*this == other);
  }

 private:
  ContemporariesIterator();

  std::unordered_set<Node*>::iterator it_set_;
  std::vector<Node*>::iterator it_vec_;
  bool use_set_;
};


class ContemporariesConstIterator {
 public:
  ContemporariesConstIterator(std::unordered_set<Node*>::const_iterator it) : it_set_(it), use_set_(true) { };
  ContemporariesConstIterator(std::vector<Node*>::const_iterator it) : it_vec_(it), use_set_(false) { };

  Node const* operator*() const {
    if (use_set_) return *it_set_;
    else return *it_vec_;
  }

  ContemporariesConstIterator& operator++() {
    if (use_set_) ++it_set_;
    else ++it_vec_;
    return *this;
  }

  bool operator==(const ContemporariesConstIterator &other) const {
    if (use_set_ != other.use_set_) return false;
    if (use_set_) return (it_set_ == other.it_set_);
    else return (it_vec_ == other.it_vec_);
  }

  bool operator!=(const ContemporariesConstIterator &other) const {
    return !(*this == other);
  }

 private:
  std::unordered_set<Node*>::const_iterator it_set_;
  std::vector<Node*>::const_iterator it_vec_;
  bool use_set_;
};

class ContemporariesContainer {
 public:
  ContemporariesContainer();
  ContemporariesContainer(const size_t pop_number,
                          const size_t sample_number,
                          RandomGenerator *rg);
  ~ContemporariesContainer(){};

  void add(Node* node);
  void remove(Node *node);
  void replaceChildren(Node *add_node);
  void replace(Node *add_node, Node *del_node_1, Node *del_node_2 = NULL);
  void clear(const bool clear_buffer = true);
  void buffer(const double current_time);

  Node* sample(const size_t pop) const;
  size_t size(const size_t pop) const;
  double buffer_time() const { return buffer_time_; }
  bool empty() const;
  bool use_set() const { return use_set_; };

  // Create Iterators
  ContemporariesConstIterator begin(const size_t pop) const {
    if (use_set_) return ContemporariesConstIterator(contemporaries_set().at(pop).cbegin());
    else return ContemporariesConstIterator(contemporaries_vector().at(pop).cbegin());
  }
  ContemporariesConstIterator end(const size_t pop) const {
    if (use_set_) return ContemporariesConstIterator(contemporaries_set().at(pop).cend());
    else return ContemporariesConstIterator(contemporaries_vector().at(pop).cend());
  }

  // Create Iterators for the buffer
  ContemporariesConstIterator buffer_begin(const size_t pop) const {
    if (use_set_) return ContemporariesConstIterator(buffer_set().at(pop).cbegin());
    else return ContemporariesConstIterator(buffer_vector().at(pop).cbegin());
  }
  ContemporariesConstIterator buffer_end(const size_t pop) const {
    if (use_set_) return ContemporariesConstIterator(buffer_set().at(pop).cend());
    else return ContemporariesConstIterator(buffer_vector().at(pop).cend());
  }
  ContemporariesIterator buffer_begin(const size_t pop) {
    if (use_set_) return ContemporariesIterator(buffer_set().at(pop).begin());
    else return ContemporariesIterator(buffer_vector().at(pop).begin());
  }
  ContemporariesIterator buffer_end(const size_t pop) {
    if (use_set_) return ContemporariesIterator(buffer_set().at(pop).end());
    else return ContemporariesIterator(buffer_vector().at(pop).end());
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

  std::vector<std::vector<Node*> > &contemporaries_vector() {
    if (use_first_) return contemporaries_vec1_;
    else return contemporaries_vec2_;
  }
  const std::vector<std::vector<Node*> > &contemporaries_vector() const {
    if (use_first_) return contemporaries_vec1_;
    else return contemporaries_vec2_;
  }
  std::vector<std::vector<Node*> > &buffer_vector() {
    if (!use_first_) return contemporaries_vec1_;
    else return contemporaries_vec2_;
  }
  const std::vector<std::vector<Node*> > &buffer_vector() const {
    if (!use_first_) return contemporaries_vec1_;
    else return contemporaries_vec2_;
  }


  std::vector<std::unordered_set<Node*> > contemporaries_set1_, contemporaries_set2_;
  std::vector<std::vector<Node*> > contemporaries_vec1_, contemporaries_vec2_;

  bool use_first_;
  bool use_set_;
  double buffer_time_;
  RandomGenerator* rg_;
};


inline ContemporariesContainer::ContemporariesContainer() {
  contemporaries_vec1_ = std::vector<std::vector<Node*> >(1, std::vector<Node*>(100));
  contemporaries_vec2_ = std::vector<std::vector<Node*> >(1, std::vector<Node*>(100));

  rg_ = NULL;
  use_first_ = true;
  buffer_time_ = -1;
  use_set_ = false;
}


inline ContemporariesContainer::ContemporariesContainer(const size_t pop_number,
                                                        const size_t sample_number,
                                                        RandomGenerator* rg) {

  // Use vectors for the storage if the number of samples is below 750.
  // This threshold is mostly arbitrary, with simulation supporting it in
  // special situations.
  if (sample_number <= 750) {
    contemporaries_vec1_ = std::vector<std::vector<Node*> >(pop_number);
    for ( auto it : contemporaries_vec1_ ) it.reserve(sample_number + 200);
    contemporaries_vec2_ = std::vector<std::vector<Node*> >(pop_number);
    for ( auto it : contemporaries_vec2_ ) it.reserve(sample_number + 200);
    use_set_ = false;
  } else {
    size_t bucket_nr = std::ceil((sample_number + 200) * 1.4);
    contemporaries_set1_ = std::vector<std::unordered_set<Node*> >(pop_number, std::unordered_set<Node*>(bucket_nr));
    contemporaries_set2_ = std::vector<std::unordered_set<Node*> >(pop_number, std::unordered_set<Node*>(bucket_nr));
    use_set_ = true;
  }
  this->rg_ = rg;
  use_first_ = true;
  buffer_time_ = DBL_MAX;
}

inline void ContemporariesContainer::add(Node* node) {
  assert(node != NULL);
  assert(!node->is_root());
  if (use_set_) contemporaries_set().at(node->population()).insert(node);
  else contemporaries_vector().at(node->population()).push_back(node);
}

inline void ContemporariesContainer::remove(Node* node) {
  assert(node != NULL);
  if (use_set_) contemporaries_set().at(node->population()).erase(node);
  else {
    size_t pop = node->population();
    auto it = std::find(contemporaries_vector().at(pop).begin(),
                        contemporaries_vector().at(pop).end(),
                        node);
    if (it != contemporaries_vector().at(pop).end()) {
      contemporaries_vector().at(pop).erase(it);
    }
  }
}

inline void ContemporariesContainer::replaceChildren(Node *add_node) {
  replace(add_node, add_node->first_child(), add_node->second_child());
}

inline void ContemporariesContainer::replace(Node *add_node, Node *del_node_1, Node *del_node_2) {
  assert(add_node != NULL);
  if (del_node_1 != NULL) remove(del_node_1);
  if (del_node_2 != NULL) remove(del_node_2);
  if (!add_node->is_root()) add(add_node);
}


inline void ContemporariesContainer::clear(const bool clear_buffer) {
  if (use_set_) {
    for (auto & pop_container : contemporaries_set()) pop_container.clear();
    if (clear_buffer) {
      for (auto & pop_container : buffer_set()) pop_container.clear();
    }
  } else {
    for (auto & pop_container : contemporaries_vector()) pop_container.clear(); 
    if (clear_buffer) {
      for (auto & pop_container : buffer_vector()) pop_container.clear();  
    }
  }
  if (clear_buffer) buffer_time_ = DBL_MAX;
  assert(this->empty());
}


inline size_t ContemporariesContainer::size(const size_t pop) const {
  if (use_set_) return contemporaries_set().at(pop).size();
  else return contemporaries_vector().at(pop).size();
}

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
  this->clear(false);
}

// Uniformly samples a random node from the current contemporaries.
// Distribution checked.
inline Node* ContemporariesContainer::sample(const size_t pop) const {
  assert( (!use_set_) || pop < contemporaries_set().size() );
  assert( (use_set_)  || pop < contemporaries_vector().size() );
  assert( this->size(pop) > 0 );

  size_t sample = rg_->sampleInt(this->size(pop));

  if (use_set_) {
    // Sample the position of the Node we return

    for (auto it = contemporaries_set().at(pop).begin();
         it != contemporaries_set().at(pop).end(); ++it) {
      assert( *it != NULL );
      if ( sample == 0 ) return (*it);
      --sample;
    }
  } else {
    return contemporaries_vector().at(pop).at(sample);
  }

  throw std::logic_error("Failed to find the contemporary I wanted to sample.");
  return NULL;
}

inline bool ContemporariesContainer::empty() const {
  if (use_set_) {
    for (auto pop_container : contemporaries_set()) {
      if ( !pop_container.empty() ) return false;
    }
  } else {
    for (auto pop_container : contemporaries_vector()) {
      if ( !pop_container.empty() ) return false;
    }
  }
  return true;
}
#endif

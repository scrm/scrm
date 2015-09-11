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
 * GNU General Public License for more details.  * * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef scrm_src_time_interval
#define scrm_src_time_interval

#include "macros.h" // Needs to be before cassert

#include "forest.h"
#include "contemporaries_container.h"
#include "model.h"

#ifdef UNITTEST
#define TI_DEBUG
#endif

#ifndef NDEBUG
#ifndef TI_DEBUG
#define TI_DEBUG
#endif
#endif

class Forest;
class TimeIntervalIterator;

class TimeInterval {
 public:
  friend class TimeIntervalIterator;

#ifdef UNITTEST
  friend class TestTimeInterval;
  friend class TestForestState;
  friend class Forest;
#endif

  TimeInterval();
  TimeInterval(TimeIntervalIterator* tii, double start_height, double end_height);
  ~TimeInterval() { };

  double start_height()  const { return this->start_height_; };
  double end_height()    const { return this->end_height_; };
  double length()        const { return (end_height() - start_height()); };
  const Forest &forest() const { return *(this->forest_); }

  //std::unordered_set<Node*> &contemporaries(const size_t pop = 0);

 private:
  double start_height_;
  double end_height_; 
  Forest* forest_;
  TimeIntervalIterator const* tii_;
};


/* WARNING: DON'T USE MULTIPLE OF THESE AT THE SAME TIME */
class TimeIntervalIterator {
 public:
  TimeIntervalIterator(Forest* forest, Node* start_node);

  void next();
  bool good() const { return this->good_; }

  TimeInterval operator*() const { return current_interval_; }
  TimeInterval operator++() { next(); return current_interval_; }
  TimeInterval operator++(int) { TimeInterval tmp = current_interval_;
                          next();
                          return tmp; }

  // Splits the current interval in two parts by adding a node inside the interval;
  // Only affects the event after the next "next()" which than represents the
  // second part of the interval.
  void splitCurrentInterval(Node* splitting_node, Node* del_node = NULL) {
    this->inside_node_ = splitting_node;
    if (del_node != NULL) contemporaries_->remove(del_node);
  };

  void recalculateInterval();

  void searchContemporaries(Node* node);
  void searchContemporariesBottomUp(Node* node, const bool use_buffer = false);

  const Forest &forest() const { return *forest_; };

#ifdef UNITTEST
  friend class TestTimeInterval;
#endif

 private:
  TimeIntervalIterator(Forest *forest);
  TimeIntervalIterator( TimeIntervalIterator const &other );
  TimeIntervalIterator& operator= ( TimeIntervalIterator const &other ); 

  Forest* forest() { return forest_; }
  ContemporariesContainer* contemporaries() { return contemporaries_; }
  Model* model() { return model_; }

  Forest* forest_;
  ContemporariesContainer* contemporaries_;
  Model* model_;

  TimeInterval current_interval_;
  double current_time_;
  NodeIterator node_iterator_;

  bool good_;
  bool model_changed_;

  Node* inside_node_;

  // Temporary values
  Node *tmp_child_1_, *tmp_child_2_, *tmp_prev_node_;
};

/** 
 * Finds all nodes which code for branches at the height of a given node in the
 * tree (i.e. the node's contemporaries). Saves this nodes in the contemporaries_
 * member.
 * 
 * @param node Node for which we are searching contemporaries
 * @return Nothing. Nodes are saved in contemporaries_.
 */
inline void TimeIntervalIterator::searchContemporaries(Node *node) {
  // Prefer top-down search if the target is above .8 coalescence units in
  // sample space. This is relatively high, but the iterative bottom-up approach
  // is faster than the recursion.
  if (node->height() >= contemporaries()->buffer_time()) {
    searchContemporariesBottomUp(node, true);
  } else {
    searchContemporariesBottomUp(node);
  }
}

#endif

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

#ifndef scrm_src_time_interval
#define scrm_src_time_interval

#include "forest.h"

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
  TimeInterval(TimeIntervalIterator const* tii, double start_height, double end_height);
  ~TimeInterval() { };

  double start_height() const { return this->start_height_; };
  double end_height()   const { return this->end_height_; };
  double length()       const { return (end_height() - start_height()); };
  const Forest &forest() const;

  size_t numberOfContemporaries(size_t pop = 0) const;
  Node* getRandomContemporary(size_t pop = 0) const;
  const std::vector<Node*> &contemporaries() const;

 private:
  Node* getIthContemporaryOfPop(size_t i, const size_t &pop) const;
  
  double start_height_;
  double end_height_; 
  TimeIntervalIterator const* tii_;
};


/* WARNING: DON'T USE MULTIPLE OF THESE AT THE SAME TIME */
class TimeIntervalIterator {
 public:
  TimeIntervalIterator();
  TimeIntervalIterator(Forest *forest, Node *start_node, bool pruning = true);
  ~TimeIntervalIterator();

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
    if (del_node != NULL) removeFromContemporaries(del_node);
  };

  void recalculateInterval();
  void removeFromContemporaries(Node* node);
  void addToContemporaries(Node* node) { 
    contemporaries_.push_back(node);
    pop_counts_.at(node->population()) += 1; 
  };
  void searchContemporariesOfNode(Node *node);
  void searchContemporariesOfNodeTopDown(Node *node, Node *current_node = NULL);
  void updateContemporaries(Node *current_node); 

  const std::vector<Node*> &contemporaries() const { return contemporaries_; };
  const Forest &forest() const { return *forest_; };
  size_t numberOfContemporaries(size_t pop = 0) const { return pop_counts_.at(pop); };

#ifdef UNITTEST
  friend class TestTimeInterval;
#endif

 private:
  TimeIntervalIterator( TimeIntervalIterator const &other );
  TimeIntervalIterator& operator= ( TimeIntervalIterator const &other ); 

  Forest* forest_;
  TimeInterval current_interval_;
  double current_time_;
  std::vector<Node*> contemporaries_;
  std::vector<size_t> pop_counts_;
  NodeIterator node_iterator_;

  bool good_;
  bool pruning_;
  bool model_changed_;

  Node* inside_node_;
};

inline const Forest &TimeInterval::forest() const { return tii_->forest(); }
#endif

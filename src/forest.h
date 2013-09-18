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

/*
 * forest.h
 * Contains the class Forest, which is the central class in scrm
 *
 * The central data structure of scrm is a forest, which is a collection of
 * trees. This class on the one hand contains a NodeContainer object with all
 * nodes building the trees, and on the other hand functions to manipulate to
 * forest.
 *
 * Most functions are defined in forest.cc with exception of pure debugging
 * functions, with are in forest-debug.cc.
 */

#ifndef scrm_src_forest
#define scrm_src_forest

//Unless compiled with options NDEBUG, we will produce a debug output using 
//'dout' instead of cout and execute (expensive) assert statements.
#ifndef NDEBUG
#define dout std::cout
#else
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cfloat>
#include <cassert>
#include <boost/assign/std/vector.hpp>
#include <ostream>

#include "node.h"
#include "event.h"
#include "model.h"
#include "param.h"
#include "node_container.h"
#include "time_interval.h"
#include "tree_point.h"
#include "random/random_generator.h"
#include "random/constant_generator.h"
#include "random/mersenne_twister.h"

class TimeInterval;
class TimeIntervalIterator;

class Forest
{
 public:

#ifdef UNITTEST
  friend class TestForest;
  friend class TestNode;
  friend class TestTimeInterval;
  friend class TestModel;
  friend class TestNodeContainer;
#endif

  friend class TimeInterval;
  friend class TimeIntervalIterator;

  Forest();
  Forest(Model *model, RandomGenerator *random_generator);
  Forest(Forest * current_forest, bool entire_ARG=true);
  virtual ~Forest();

  //Getters & Setters
  const Model &model() const { return *model_; }
  void set_model(Model* model) { this->model_ = model; }
  Model* writable_model() { return this->model_; };

  Node* local_root() const { return local_root_; }
  void set_local_root(Node* local_root) { local_root_ = local_root; };

  Node* primary_root() const { return primary_root_; }
  void set_primary_root(Node* primary_root) { primary_root_ = primary_root; };

  size_t sample_size() const { 
    if (sample_size_ == 0) return model().sample_size(); 
    return this->sample_size_; 
  }
  void set_sample_size(const size_t &size ) { sample_size_ = size; }

  double current_base() const { return current_base_; }
  void set_current_base(double base) { current_base_ = base; }

  double next_base() const{return next_base_;}
  void set_next_base() {
    next_base_ = current_base_ + random_generator()->sampleExpo(local_tree_length() * model().recombination_rate());
  } 

  double local_tree_length() const { return this->local_root()->length_below(); }

  void set_random_generator(RandomGenerator *rg) {
    this->random_generator_ = rg; }
  RandomGenerator* random_generator() const { return this->random_generator_; }

  NodeContainer const *getNodes() const { return nodes_; };

  size_t prune_countdown() const{return prune_countdown_;}  // We will prune once this countdown reaches 0
  void set_prune_countdown(size_t  prune_countdown){prune_countdown_=prune_countdown_;}

  bool pruning() const{return pruning_;}
  void set_pruning(bool pruning){pruning_=pruning;}

  // Central functions
  void buildInitialTree();
  void sampleNextGenealogy();

  //Debugging Tools
  void addNodeToTree(Node *node, Node *parent, Node *first_child, Node *second_child);
  void createExampleTree();
  bool checkLeafsOnLocalTree(Node const* node=NULL) const;
  bool checkTree(Node const* root = NULL) const;
  double calcTreeLength() const;
  bool checkTreeLength() const;
  bool checkInvariants(Node const* node = NULL) const;
  bool checkNodeProperties() const;
  bool printNodes() const;

  //Tree printing
  int countLinesLeft(Node const* node) const;
  int countLinesRight(Node const* node) const;
  int countBelowLinesLeft(Node const* node) const;
  int countBelowLinesRight(Node const* node) const;
  bool printTree();
  bool printTree_cout();

  std::vector<Node const*> determinePositions() const;
  void printPositions(const std::vector<Node const*> &positions) const;

  NodeContainer *nodes() { return this->nodes_; }

  //segegrating sites
  std::ostream &generateSegData(std::ostream &output, int total_mut);

  TreePoint samplePoint(Node* node = NULL, double length_left = -1);

  //derived class from Forest

  //virtual void initialize_event(const TimeInterval & current_event, double current_time);
  virtual	void initialize_recomb_coalescent(const double rec_height);
  //virtual	void record_coalevent();
  //virtual	void record_recombevent();
  virtual void initialize_event(double start_time);
  virtual void record_coalevent(const TimeInterval & current_event, double end_time);
  virtual void record_recombevent(const TimeInterval & current_event, double end_time);
  virtual	void clear_initial_coalevent();

 private:
  //segegrating sites  
  void find_descndnt();
  void exp_mut_num(int total_mut);

  //Operations on the Tree
  Node* cut(const TreePoint &cut_point);
  void updateAbove(Node* node, 
                   bool above_local_root = false,
                   bool recursive = true,
                   bool local_only = false);

  // Tools for doing coalescence & recombination
  void sampleCoalescences(Node *start_node, bool pruning);
  //making samplePoint public
  //TreePoint samplePoint(Node* node = NULL, double length_left = -1);
  size_t getNodeState(Node const *node, const double &current_time) const;
  Node* updateBranchBelowEvent(Node* node, const TreePoint &event_point); 

  Node* possiblyMoveUpwards(Node* node, const TimeInterval &event);

  // Implementation of the different events
  void implementCoalescence(const Event &event, TimeIntervalIterator &tii);
  void implementPwCoalescence(Node* root_1, Node* root_2, const double &time);
  void implementRecombination(const Event &event, TimeIntervalIterator &tii);
  void implementMigration(const Event &event, TimeIntervalIterator &tii);
  void implementFixedTimeEvent(TimeIntervalIterator &ti);

  // Pruning
  bool isPrunable(Node const* node) const;
  void prune(Node* node); 

  // Calculation of Rates
  double calcCoalescenceRate(const size_t &pop, const TimeInterval &ti) const;

  double calcPwCoalescenceRate(const size_t &pop) const {
    // Rate a pair is 1/(2N), as N is the diploid population size
    return ( 1.0 / ( 2.0 * this->model().population_size(pop) ) );
  }

  double calcRecombinationRate(Node const* node) const {
    return ( model().recombination_rate() * (this->current_base() - node->last_update()) );
  }

  void calcRates(const TimeInterval &ti);

  void sampleEvent(const TimeInterval &ti, double tmp_event_time,  
                   size_t tmp_event_line, Event &return_event) const;

  void sampleEventType(const double &time, const size_t &time_line, 
                       const TimeInterval &ti, Event &return_event) const;

  void selectFirstTime(const double &new_time, const size_t &time_line, 
                       double &current_time, size_t &current_time_line) const;

  double getTimeLineGrowth(const size_t &time_line) const {
    if (time_line == 0) return 0.0;
    else if (time_line == 1) return model().growth_rate(active_node(0)->population());
    else if (time_line == 2) return model().growth_rate(active_node(1)->population());
    else throw std::out_of_range("Trying to get growthrate of unknown time line.");
  }


  //Private variables
  NodeContainer* nodes_;    // The nodes of the Tree/Forest

  // We have 2 different roots that are important:
  // - local_root: root of the smallest subtree containing all local sequences
  // - primary_root: root of the tree that contains all local sequences (do we
  //   need this one?)

  Node* local_root_;
  Node* primary_root_;

  double current_base_;     // The current position of the sequence we are simulating
  double next_base_;
  size_t sample_size_;      // The number of sampled nodes (changes while building the initial tree)
  size_t prune_countdown_;  // We will prune once this countdown reaches 0
  bool pruning_;

  Model* model_;
  RandomGenerator* random_generator_;

  void initialize(Model *model = new Model(),
                  RandomGenerator *rg = NULL);

  void createSampleNodes();

  // Temporarily used when doing the coalescence

  // Rates: 
  double rates_[3];

  // States: Each (branch above an) active node can either be in state
  // - 0 = off (the other coalescence has not reached it yet) or
  // - 1 = potentially coalescing in a time interval or
  // - 2 = potentially recombining in a time interval
  size_t states_[2];
  Node* active_nodes_[2];
  Event  tmp_event_;
  size_t tmp_event_line_;
  double tmp_event_time_;

  // These are pointers to the up to two active nodes during a coalescence
  size_t active_nodes_timelines_[2];
  Node* active_node(size_t nr) const { return active_nodes_[nr]; };
  void set_active_node(size_t nr, Node* node) { active_nodes_[nr] = node; };

  Node* getEventNode() const { return active_node(tmp_event_.active_node_nr()); }
  size_t getEventNodesState() const { return states_[tmp_event_.active_node_nr()]; };
  Node* getOtherNode() const { return active_node(1-tmp_event_.active_node_nr()); };
  size_t getOtherNodesState() const { return states_[1-tmp_event_.active_node_nr()]; };

  bool coalescence_finished_;
};

bool areSame(const double &a, const double &b, 
             const double &epsilon = std::numeric_limits<double>::epsilon());

#endif

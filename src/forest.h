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

/*!
 * \file forest.h
 * \brief Contains the class Forest, which is the central class in scrm
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

#include <string>
#include <vector>
#include <unordered_set>
#include <stdexcept>
#include <cassert>
#include <iostream> // ostreams
#include <iomanip>  // Used for debug output
#include <sstream>  // Used for debug output

#include "contemporaries_container.h"
#include "event.h"
#include "model.h"
#include "macros.h"
#include "node.h"
#include "node_container.h"
#include "time_interval.h"
#include "tree_point.h"
#include "random/random_generator.h"
#include "summary_statistics/summary_statistic.h"
extern int recombination_counter; //DEBUG

class TimeInterval;
class TimeIntervalIterator;
enum eventCode { NOEVENT, EVENT, INIT_NULL};

class Forest
{
 public:
  Node* readNewickNode( std::string &in_str, std::string::iterator &current_it, size_t parenthesis_balance = 0, Node* parent = NULL );
  void readNewick(std::string &in_str);
  ContemporariesContainer* contemporaries()  {return &this->contemporaries_;};

#ifdef UNITTEST
  friend class TestForest;
  friend class TestNode;
  friend class TestTimeInterval;
  friend class TestModel;
  friend class TestNodeContainer;
#endif
  friend class ForestState;
  friend class ParticleContainer;
  friend class TimeInterval;
  friend class TimeIntervalIterator;

  Forest(Model *model, RandomGenerator *random_generator);
  Forest(const Forest &current_forest);
  virtual ~Forest() {};

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
  void set_sample_size(const size_t size ) { sample_size_ = size; }

  size_t segment_count() const { return current_rec_; }

  void sampleNextBase();

  // Resamples next recombination base, in case a state was duplicated
  // Note that current_base_ must be up to date
  void resampleNextBase() {
      // first undo any (now possibly premature) change in the recombination rate
      if (next_base() == model().getCurrentSequencePosition()) {
          writable_model()->decreaseSequencePosition();
      }
      // the remainder is the same as for sampleNextBase; the previous next_base_ is overwritten
      sampleNextBase();
  }

  /**
   * @brief Returns the length of the sequence for with the current tree is
   * valid
   *
   * @param finite_sites If 'true', the length is measured in number of bases
   * (e.g. integer sequence positions) for which the tree is valid. Otherwise,
   * the length is measured real valued number on a continuous chromosome.  
   *
   * @return The length of the current segment (see above for its unit)
   */
  double calcSegmentLength() const {
    if (model().getSequenceScaling() == relative) {
      return (next_base() - current_base()) / model().loci_length();
    } else {
      return ceil(next_base()) - ceil(current_base());
    }
  }

  void set_random_generator(RandomGenerator *rg) {
    this->random_generator_ = rg; 
  }
  RandomGenerator* random_generator() const { return this->random_generator_; }

  NodeContainer const *getNodes() const { return &nodes_; };

  // Central functions
  void buildInitialTree();
  void sampleNextGenealogy( bool recordEvents = false );
  TreePoint samplePoint(Node* node = NULL, double length_left = -1) const;

  void clear();

  //Debugging Tools
  void addNodeToTree(Node *node, Node *parent, Node *first_child, Node *second_child);
  void createExampleTree();
  void createScaledExampleTree();
  bool checkLeafsOnLocalTree(Node const* node=NULL) const;
  bool checkTree(Node const* root = NULL) const;
  double calcTreeLength() const;
  bool checkTreeLength() const;
  bool checkInvariants(Node const* node = NULL) const;
  bool checkNodeProperties() const;
  bool checkContemporaries(const double time) const;
  bool printNodes() const;
  bool checkForNodeAtHeight(const double height) const;
  bool checkRootIsRegistered(Node const* node) const;
  bool checkRoots() const;

  //Debug Tree Printing
  int countLinesLeft(Node const* node) const;
  int countLinesRight(Node const* node) const;
  int countBelowLinesLeft(Node const* node) const;
  int countBelowLinesRight(Node const* node) const;
  bool printTree() const;
  void printTree_cout() const;
  std::vector<Node const*> determinePositions() const;
  void printPositions(const std::vector<Node const*> &positions) const;

  NodeContainer *nodes() { return &(this->nodes_); }

  double getTMRCA(const bool &scaled = false) const {
    if (scaled) return local_root()->height() / (4 * this->model_->default_pop_size);
    else return local_root()->height();
  }

  double getLocalTreeLength(const bool &scaled = false) const {
    if (scaled) return local_root()->length_below() / (4 * this->model_->default_pop_size);
    else return local_root()->length_below();
  }

  Node* trackLocalNode(Node *node) const;

  //record events
  virtual void record_Recombevent_b4_extension ( ) { }
  virtual void record_Recombevent_atNewGenealogy ( double event_height ){ }
  virtual void record_all_event( TimeInterval const &ti, double &recomb_opp_x_within_scrm ) { }
  double recomb_opp_x_within_scrm; // DEBUG

  // Calc & Print Summary Statistics
  void calcSegmentSumStats();
  void clearSumStats();
  void printLocusSumStats(std::ostream &output) const;
  void printSegmentSumStats(std::ostream &output) const;

  double get_rec_base(const size_t idx) const {
    assert(idx < rec_bases_.size());
    return rec_bases_[idx];
  }

  double current_base() const { return get_rec_base(current_rec_); }
  double next_base() const { return get_rec_base(current_rec_ + 1); }
  void set_current_base(double const base) { rec_bases_[current_rec_] = base; }; 
  void set_next_base(double const base) { rec_bases_.push_back(base); }; 

  size_t current_rec() const { return current_rec_; };

 private:
  Forest() { this->initialize(); }

  //Operations on the Tree
  Node* cut(const TreePoint &cut_point);

  void updateAbove(Node* node, 
                   bool above_local_root = false,
                   const bool &recursive = true,
                   const bool &invariants_only = false);

  // Tools for doing coalescence & recombination
  void sampleCoalescences(Node *start_node,  bool recordEvents = true );
  size_t getNodeState(Node const *node, const double current_time) const;
  Node* updateBranchBelowEvent(Node* node, const TreePoint &event_point); 
  Node* possiblyMoveUpwards(Node* node, const TimeInterval &event);

  // Implementation of the different events
  void implementNoEvent(const TimeInterval &ti, bool &coalescence_finished);
  void implementCoalescence(const Event &event, TimeIntervalIterator &tii);
  void implementPwCoalescence(Node* root_1, Node* root_2, const double time);
  void implementRecombination(const Event &event, TimeIntervalIterator &tii);
  void implementMigration(const Event &event, const bool &recalculate, TimeIntervalIterator &tii);
  void implementFixedTimeEvent(TimeIntervalIterator &ti);

  // Pruning
  bool nodeIsOld(Node const* node) const {
    if ( node->local() ) return false;
    if ( node->is_root() ) return false;
    if ( model().has_window_rec() && 
         segment_count() - node->last_update() > model().window_length_rec()) {
      return true;
    }
    if ( model().has_window_seq() && 
         current_base() - get_rec_base(node->last_update()) > model().window_length_seq()) {
      return true;
    }
    return false;
  }

  bool nodeIsActive(Node const* node) const {
    return (node == active_node(0) || node == active_node(1));
  }

  bool pruneNodeIfNeeded(Node* node, const bool prune_orphans = true);

  // Calculation of Rates
  double calcCoalescenceRate(const size_t pop, const TimeInterval &ti) const;
  double calcPwCoalescenceRate(const size_t pop, const TimeInterval &ti) const;
  double calcRecombinationRate(Node const* node) const;
  void calcRates(const TimeInterval &ti, double &recomb_opp_x_within_scrm);

  void sampleEvent(const TimeInterval &ti, double &event_time, Event &return_event) const;

  void sampleEventType(const double time, const size_t time_line, 
                       const TimeInterval &ti, Event &return_event) const;

  void selectFirstTime(const double new_time, const size_t time_line, 
                       double &current_time, size_t &current_time_line) const;

  double getTimeLineGrowth(const size_t time_line) const {
    if (time_line == 0) return 0.0;
    else if (time_line == 1) return model().growth_rate(active_node(0)->population());
    else if (time_line == 2) return model().growth_rate(active_node(1)->population());
    else throw std::out_of_range("Trying to get growthrate of unknown time line.");
  }

  // Private Members
  NodeContainer nodes_;    // The nodes of the Tree/Forest

  // We have 3 different kind of roots that are important:
  // local root: root of the smallest subtree containing all local sequences
  Node* local_root_;

  // primary root: root of the tree that contains all local sequences
  Node* primary_root_;

  // secondary roots: roots of trees that contain only non-local nodes
  // std::unordered_set<Node*> secondary_roots_;

  size_t sample_size_;      // The number of sampled nodes (changes while building the initial tree)
  size_t current_rec_;      // A counter for recombinations

  std::vector<double> rec_bases_; // Genetic positions of the recombinations 

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
  double tmp_event_time_;
  ContemporariesContainer contemporaries_;

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

#endif

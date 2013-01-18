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
#include <iostream>
#include <stdexcept>
#include <cfloat>
#include <cassert>
#include <boost/assign/std/vector.hpp>

#include "node.h"
#include "model.h"
#include "node_container.h"
#include "event.h"
#include "tree_point.h"
#include "random/random_generator.h"
#include "random/constant_generator.h"
#include "random/fake_generator.h"
#include "random/mersenne_twister.h"

class Forest
{
 public:

#ifdef UNITTEST
  friend class TestForest;
  friend class TestNode;
#endif
  friend class Event;
  friend class EventIterator;

  Forest();
  Forest(Model model, RandomGenerator *random_generator);
  ~Forest();

  //Getters & Setters
  Model model() const { return this->model_; }
  void set_model(const Model &model) { this->model_ = model; }

  Node* ultimate_root() const { return ultimate_root_; }
  Node* local_root() const { return local_root_; }
  void set_local_root(Node* local_root) { local_root_ = local_root; };
  
  Node* primary_root() const { return primary_root_; }
  void set_primary_root(Node* primary_root) { primary_root_ = primary_root; };

  size_t sample_size() const { return this->sample_size_; }

  size_t current_base() const { return current_base_; }
  void set_current_base(size_t base) { current_base_ = base; }

  double local_tree_length() const { return this->local_root()->length_below(); }
  RandomGenerator* random_generator() const { return this->random_generator_; }
  NodeContainer const *getNodes() const { return nodes_; };

  //Operations on the Tree
  void addNodeToTree(Node *node, Node *parent, Node *lower_child, Node *higher_child);
  void cut(const TreePoint &cut_point);
  void deactivateSubtree(Node* node);
  void updateAbove(Node* node, 
                   bool above_local_root = false,
                   bool recursive = true,
                   bool local_only = false);

  Node* moveUpwardsInTree(Node* node);


  //Operations to manage the fake binary tree connecting all trees of the forest
  void createRoots();
  void registerNonLocalRoot(Node *node, Node *position = NULL);
  void unregisterNonLocalRoot(Node* node, Node *position = NULL);

  //
  void buildInitialTree();
  TreePoint samplePoint(Node* node = NULL, double length_left = -1);
  void sampleNextGenealogy();
  void sampleCoalescences(Node *start_node, const bool &for_initial_tree = false);
  double calcCoalescenceRate(int lines_number, int coal_lines_number = 1);
  void coalesNodeIntoTree(Node* coal_node, const TreePoint &coal_point);

  //Debugging Tools
  void createExampleTree();
  bool checkLeafsOnLocalTree(Node const* node=NULL) const;
  bool checkTree(Node* root = NULL) const;
  double calcTreeLength() const;
  bool checkTreeLength() const;
  bool checkNodeProperties() const;
  bool printNodes() const;

  //Tree printing
  int countLinesLeft(Node const* node) const;
  int countLinesRight(Node const* node) const;
  int countBelowLinesLeft(Node const* node) const;
  int countBelowLinesRight(Node const* node) const;
  bool printTree() const;


 private:
  //Private variables
  NodeContainer* nodes_;    // The nodes of the Tree/Forest

  // We have 3 different roots that are important:
  // - local_root: root of the smallest subtree containing all local sequences
  // - primary_root: root of the tree that contains all local sequences
  // - ultimate root: root of the fake tree structure connecting all trees of
  //                the forest into a tree.

  Node* local_root_;
  Node* primary_root_;
  Node* ultimate_root_;      

  size_t current_base_;     // The current position of the sequence we are simulating

  size_t sample_size_;      // The number of sampled nodes (changes while building the initial tree)

  Model model_;
  RandomGenerator* random_generator_;

  double local_tree_length_;
  double total_tree_length_;
  
  void initialize(Model model = Model(),
                  RandomGenerator* rg = NULL, 
                  Node* ultimate_root = NULL);

  NodeContainer *nodes() { return this->nodes_; }
  NodeContainer const *nodes() const { return this->nodes_; }

  void set_random_generator(RandomGenerator *rg) {
    this->random_generator_ = rg; }

  void set_sample_size(const size_t &size ) { sample_size_ = size; }

  void set_ultimate_root(Node* ultimate_root) { ultimate_root_ = ultimate_root; }

  void createSampleNodes();
};

bool areSame(double a, double b);

#endif

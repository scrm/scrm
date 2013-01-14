#ifndef scrm_src_forest
#define scrm_src_forest

//Unless compiled with options NDEBUG, we will produce a debug output using 
//'dout' instead of cout and execute (expensive) assert statements.
#ifndef NDEBUG
#define dout std::cout
#else
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
  Node* local_root() const { return ultimate_root()->lower_child(); }
  int sample_size() const { return this->model().sample_size(); }
  double local_tree_length() const { return this->local_tree_length_; }
  double total_tree_length() const { return this->total_tree_length_; }
  RandomGenerator* random_generator() const { return this->random_generator_; }
  NodeContainer const *getNodes() const { return nodes_; };

  //Operations on the Tree
  void addNodeToTree(Node *node, Node *parent, Node *lower_child, Node *higher_child);
  void cut(const TreePoint &cut_point);

  void calcTreeLength(double *local_length, double *total_length) const;
  void updateTreeLength();

  //Operations to manage the fake binary tree connecting all trees of the forest
  void createRoots();
  void registerNonLocalRoot(Node *node, Node *position = NULL);
  void unregisterNonLocalRoot(Node* node, Node *position = NULL);

  //
  void buildInitialTree();
  void buildInitialTree_old();
  TreePoint samplePoint(bool only_local = false);
  void sampleNextGenealogy();
  void sampleCoalescences(Node *start_node, const bool &for_initial_tree = false);
  double calcCoalescenceRate(int lines_number, int coal_lines_number = 1);
  void coalesNodeIntoTree(Node* coal_node, const TreePoint &coal_point);

  //Debugging Tools
  void createExampleTree();
  bool checkTree(Node* root = NULL) const;
  bool checkTreeLength() const;
  bool checkNodesSorted() const;
  bool printNodes() const;

  //Tree printing
  int countLinesLeft(Node const* node) const;
  int countLinesRight(Node const* node) const;
  int countBelowLinesLeft(Node const* node) const;
  int countBelowLinesRight(Node const* node) const;
  bool printTree() const;


 private:
  //Private variables
  NodeContainer* nodes_;
  Node* ultimate_root_;

  Model model_;
  RandomGenerator* random_generator_;

  double local_tree_length_;
  double total_tree_length_;
  
  void initialize(Model model = Model(),
                  RandomGenerator* rg = NULL, 
                  Node* ultimate_root = NULL, 
                  int local_tree_length = 0,
                  int total_tree_length = 0);

  NodeContainer *nodes() { return this->nodes_; }
  NodeContainer const *nodes() const { return this->nodes_; }

  void set_random_generator(RandomGenerator *rg) {
    this->random_generator_ = rg; }

  void set_local_tree_length(const double &length) { local_tree_length_ = length; }
  void set_total_tree_length(const double &length) { total_tree_length_ = length; }
  void set_ultimate_root(Node* ultimate_root) { ultimate_root_ = ultimate_root; }

  void inc_tree_length(const double &by, const bool &active);
  void dec_tree_length(const double &by, const bool &active) 
    { inc_tree_length(-1 * by, active); }

  void createSampleNodes();
};

bool areSame(double a, double b);

#endif

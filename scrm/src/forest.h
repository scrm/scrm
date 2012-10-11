#ifndef scrm_src_forest
#define scrm_src_forest

#include <vector>
#include <iostream>
#include <stdexcept>
#include <cfloat>
#include "node.h"
#include "model.h"
#include "event.h"
#include "random/random_generator.h"

class Forest
{
 public:

#ifdef UNITTEST
  friend class TestForest;
#endif

  Forest();
  Forest(Model model, RandomGenerator *random_generator);
  ~Forest();

  //Getters & Setters
  Model model() { return this->model_; }
  void set_model(const Model &model) { this->model_ = model; }
  Node* ultimate_root() { return ultimate_root_; }
  void set_ultimate_root(Node* ultimate_root) { ultimate_root_ = ultimate_root; }
  int sample_size() { return this->model().sample_size(); }
  double local_tree_length() { return this->local_tree_length_; }
  void set_local_tree_length(const double &length) { local_tree_length_ = length; }
  double total_tree_length() { return this->total_tree_length_; }
  void set_total_tree_length(const double &length) { total_tree_length_ = length; }

  //Operations on Nodes
  Node* getFirstNode();
  vector<Node*>::iterator getNodesEnd();
  vector<Node*>::iterator getNodeFwIterator();

  void addNode(Node *node);
  void addNodeAfter(const Node &node, const Node &after_node);
  void addNodeBefore(const Node &node, const Node &before_node);
  int countNodes();

  //Operations on the Tree
  void addNodeToTree(Node *node, Node *parent, Node *lower_child, Node *higher_child);

  //
  void buildInitialTree();
  void samplePoint(bool only_local, Node **above_node, double *height_above);
  void sampleNextGenealogy();

  //Debugging Tools
  void createExampleTree();
  void checkTree(Node* root = NULL);
  void checkNodesSorted();
  void printNodes();

 private:
  //Private variables
  std::vector<Node*> nodes_;
  Node* ultimate_root_;
  Model model_;
  RandomGenerator* random_generator_;
  
  void initialize(Model model = Model(),
                  RandomGenerator* rg = NULL, 
                  Node* ultimate_root = NULL, 
                  int local_tree_length = 0,
                  int total_tree_length = 0);



  double local_tree_length_;
  double total_tree_length_;

  std::vector<Node*> nodes() { return this->nodes_; }

  RandomGenerator* random_generator() { return this->random_generator_; }
  void set_random_generator(RandomGenerator *rg) {
    this->random_generator_ = rg; }

  void inc_local_tree_length(const double &by);
  void dec_local_tree_length(const double &by) { inc_local_tree_length(-1 * by); }
  void inc_total_tree_length(const double &by);
  void dec_total_tree_length(const double &by) { inc_total_tree_length(-1 * by); }

  void createSampleNodes();
  
};

#endif

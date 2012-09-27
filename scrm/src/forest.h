#ifndef scrm_src_forest
#define scrm_src_forest

#include <vector>
#include <iostream>
#include "node.h"
#include "model.h"
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
   int sample_size() { return this->model().sample_size(); }
   int local_tree_length() { return this->local_tree_length_; }
   int total_tree_length() { return this->total_tree_length_; }

   //Operations on Nodes
   void addNode(Node *node);
   void addNodeAfter(const Node &node, const Node &after_node);
   void addNodeBefore(const Node &node, const Node &before_node);
   int countNodes();
   Node* getFirstNode();

   //
   void buildInitialTree();
   void samplePoint(bool only_local, Node **above_node, double *height_above);
   void addNodeToTree(Node *node, Node *parent, Node *lower_child, Node *higher_child);

  private:
   std::vector<Node*> nodes_;
   Node* ultimate_root_;
   Model model_;
   RandomGenerator *random_generator_;

   int local_tree_length_;
   int total_tree_length_;

   std::vector<Node*> nodes() { return this->nodes_; }

   RandomGenerator* random_generator() { return this->random_generator_; }
   void set_random_generator(RandomGenerator *rg) {
     this->random_generator_ = rg; }

   void inc_local_tree_length(const double &by);
   void dec_local_tree_length(const double &by) { inc_local_tree_length(-1 * by); }
   void inc_total_tree_length(const double &by);
   void dec_total_tree_length(const double &by) { inc_total_tree_length(-1 * by); }

   void createSampleNodes();
   void printNodes();
};

#endif

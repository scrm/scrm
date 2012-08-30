#ifndef scrm_src_forest
#define scrm_src_forest

#include <vector>
#include <iostream>
#include "node.h"
#include "model.h"
#include "random.h"

class Forest
{
  public:

#ifdef UNITTEST
    friend class TestForest;
#endif

   Forest();
   Forest(Model model, RandomGenerator random_generator);
   ~Forest();
                         
   //Getters & Setters
   Model model() { return this->model_; }
   void set_model(const Model &model) { this->model_ = model; }
   int sample_size() { return this->model().sample_size(); }

   //Operations on Nodes
   void addNode(Node *node);
   void addNodeAfter(const Node &node, const Node &after_node);
   void addNodeBefore(const Node &node, const Node &before_node);

   int countNodes();

   //
   void buildInitialTree();

  private:
   std::vector<Node*> nodes_;
   Model model_;
   RandomGenerator random_generator_;

   std::vector<Node*> nodes() { return this->nodes_; }

   RandomGenerator* random_generator() { return &(this->random_generator_); }
   void set_random_generator(RandomGenerator rg) {
     this->random_generator_ = rg; }

   void createSampleNodes();
   void printNodes();
};

#endif

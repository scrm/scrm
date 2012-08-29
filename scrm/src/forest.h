#ifndef scrm_src_forest
#define scrm_src_forest

#include <vector>
#include "node.h"
#include "model.h"

class Forest
{
  public:

#ifdef UNITTEST
    friend class TestForest;
#endif

   Forest();
   Forest(Model model);
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

   std::vector<Node*> nodes() { return this->nodes_; }

   void createSampleNodes();
   void printNodes();
};

#endif

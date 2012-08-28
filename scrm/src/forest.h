#ifndef scrm_forest
#define scrm_forest

#include <vector>
#include "node.h"

class Forest
{
  public:                       
   Forest();                      
   ~Forest();
                         
   //Getters & Setters

   //Operations on Nodes
   void addNode(const Node &node);
   void addNodeAfter(const Node &node, const Node &after_node);
   void addNodeBefore(Node *node, Node *before_node);

   int countNodes();

  private:
   std::vector<Node> nodes_;
};

#endif

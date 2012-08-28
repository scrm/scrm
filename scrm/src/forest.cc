#include "forest.h"

Forest::Forest() { };
Forest::~Forest() { };

void Forest::addNode(const Node &node) {
  nodes_.push_back(node);
}

void addNodeAfter(const Node &node, const Node &after_node){
  
}

int Forest::countNodes(){
  return(nodes_.size());
}

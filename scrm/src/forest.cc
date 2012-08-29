#include "forest.h"
#include "random.h"
#include <vector>
#include <iostream>

using namespace std;

Forest::Forest() {

};

Forest::Forest(Model model) {
  this->set_model(model);
}

Forest::~Forest() { };

void Forest::addNode(Node *node) {
  nodes_.push_back(node);
}

void Forest::addNodeAfter(const Node &node, const Node &after_node){
  
}

void Forest::addNodeBefore(const Node &node, const Node &before_node){

}

//vector<Node*> Forest::nodes() {
//  vector<Node*> nodes;
//  nodes.reserve(this->nodes_.size());
//  for(vector<Node>::iterator it = this->nodes_.begin(); it != this->nodes_.end(); ++it) {
//        nodes.push_back(&(*it));
//  }
//  return(nodes);
//}

void Forest::printNodes() {
  for(int i = 0; i < this->countNodes(); ++i) {
        cout << "Addr: " << this->nodes()[i] << " | ";
        cout << "Height: " << this->nodes()[i]->height() << " | ";
        cout << "Parent: " << this->nodes()[i]->parent() << endl;
  }
}


void Forest::buildInitialTree() {
  this->createSampleNodes();

  RandomGenerator rg = RandomGenerator();

  cout << "Preparing coalescence" << endl;
  vector<Node*> uncoalesced_nodes = this->nodes();
  double time = 0;
  while (uncoalesced_nodes.size() > 1) {
    int node1, node2;
    int n = uncoalesced_nodes.size();
    double rate = (1.0/(2*this->model().population_size()))*n*(n-1)/2.0;
    time += rg.sampleExpo(rate);
    rg.sampleTwoElements(uncoalesced_nodes.size(), &node1, &node2);
    cout << "Coalescing Nodes " << node1 << " and " << node2 << endl;

    //Creating parent
    Node* parent = new Node(time);
    uncoalesced_nodes.push_back(parent);
    this->addNode(parent);

    //Implementing
    uncoalesced_nodes[node1]->set_parent(parent);
    uncoalesced_nodes[node2]->set_parent(parent);

    //Erase coalesced nodes
    uncoalesced_nodes.erase(uncoalesced_nodes.begin() + node1);
    if (node1 < node2) --node2;
    uncoalesced_nodes.erase(uncoalesced_nodes.begin() + node2);

    this->printNodes();
  }
  cout << uncoalesced_nodes.size() << endl;
}

void Forest::createSampleNodes() {
  for (int i=0; i < this->sample_size(); i++) {
    Node* node = new Node(0);
    this->addNode(node);
  }
}

int Forest::countNodes(){
  return(nodes_.size());
}


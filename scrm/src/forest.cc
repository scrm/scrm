#include "forest.h"

using namespace std;

Forest::Forest() {

};

Forest::Forest(Model model, RandomGenerator *random_generator) {
  this->set_model(model);
  this->set_random_generator(random_generator);
}

Forest::~Forest() { };

void Forest::addNode(Node *node) {
  nodes_.push_back(node);
}

void Forest::addNodeAfter(const Node &node, const Node &after_node){
  
}

void Forest::addNodeBefore(const Node &node, const Node &before_node){

}

void Forest::addNodeToTree(Node *node, Node *parent, Node *lower_child, Node *higher_child) {
  this->addNode(node);
  
  if (parent != NULL) {
    node->set_parent(parent);
    if (parent->lower_child() == NULL) parent->set_lower_child(node);
    else {
      if (parent->lower_child()->height() > node->height()) {
        parent->set_higher_child(parent->lower_child());
        parent->set_lower_child(node);
      } else {
        parent->set_higher_child(node);
      }
    }
    double height_above = parent->height() - node->height();
    node->set_height_above(height_above);
    this->inc_local_tree_length(height_above);
    this->inc_total_tree_length(height_above);
  }

  if (lower_child != NULL) {
    node->set_lower_child(lower_child);
    lower_child->set_parent(node);
    
    double height_above = node->height() - lower_child->height();
    lower_child->set_height_above(height_above);
    this->inc_local_tree_length(height_above);
    this->inc_total_tree_length(height_above);
  }

  if (higher_child != NULL) {
    node->set_higher_child(higher_child);
    higher_child->set_parent(node);

    double height_above = node->height() - higher_child->height();
    higher_child->set_height_above(height_above);
    this->inc_local_tree_length(height_above);
    this->inc_total_tree_length(height_above);
  }
}

void Forest::printNodes() {
  for(int i = 0; i < this->countNodes(); ++i) {
        cout << "Addr: " << this->nodes()[i] << " | ";
        cout << "Height: " << this->nodes()[i]->height() << " | ";
        cout << "Parent: " << this->nodes()[i]->parent() << endl;
  }
}

Node* Forest::getFirstNode() {
  return(this->nodes_[0]);
}

void Forest::buildInitialTree() {
  this->createSampleNodes();

  RandomGenerator *rg = this->random_generator();
  cout << "Preparing coalescence" << endl;
  vector<Node*> uncoalesced_nodes = this->nodes();
  double time = 0;
  while (uncoalesced_nodes.size() > 1) {
    int node1, node2;
    int n = uncoalesced_nodes.size();
    double rate = (1.0/(2*this->model().population_size()))*n*(n-1)/2.0;
    time += rg->sampleExpo(rate);
    rg->sampleTwoElements(uncoalesced_nodes.size(), &node1, &node2);
    cout << "Coalescing Nodes " << node1 << " and " << node2 << endl;

    //Creating parent and add it to the tree
    Node* parent = new Node(time);
    uncoalesced_nodes.push_back(parent);
    this->addNodeToTree(parent, NULL, uncoalesced_nodes[node1], uncoalesced_nodes[node2]);

    //Erase coalesced nodes
    uncoalesced_nodes.erase(uncoalesced_nodes.begin() + node1);
    if (node1 < node2) --node2;
    uncoalesced_nodes.erase(uncoalesced_nodes.begin() + node2);

    this->printNodes();
  }
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

void Forest::inc_local_tree_length(const double &by) {
  local_tree_length_ = local_tree_length_ + by;
}

void Forest::inc_total_tree_length(const double &by) {
  total_tree_length_ = total_tree_length_ + by;
}


void Forest::samplePoint(bool only_local = true, Node** above_node=NULL, double* height_above=NULL) {
  //O(#Nodes) implementation
  double length = 0;

  if (only_local) length = this->local_tree_length();
  else            length = this->total_tree_length();

  double point = this->random_generator()->sample() * length;
  cout << "Sampled at height " << point << endl;
  
  for (vector<Node*>::iterator it = nodes_.begin(); it!=nodes_.end(); ++it) {
        cout << (*it)->height_above() << endl;
        if ((*it)->height_above() > point) {
          cout << "HERE" << endl;
          above_node = &(*it);
          *height_above = point;
          break;
        }
        point -= (*it)->height_above();
        cout << "LEFT: " << point << endl;
  }
}

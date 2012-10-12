#include "forest.h"


/******************************************************************
 * Constructors & initialization
 *****************************************************************/

Forest::Forest() {
  this->initialize();
};

Forest::Forest(Model model, RandomGenerator* random_generator) {
  this->initialize(model, random_generator);
}

// Sets member variable to default values
void Forest::initialize(Model model, 
                        RandomGenerator* rg, 
                        Node* ultimate_root, 
                        int local_tree_length,
                        int total_tree_length) {

  this->set_model(model);
  this->set_random_generator(rg);
  this->set_ultimate_root(ultimate_root);
  this->set_local_tree_length(local_tree_length);
  this->set_total_tree_length(total_tree_length);
}



/******************************************************************
 * Destructor
 *****************************************************************/

Forest::~Forest() { 
  for (std::vector<Node*>::iterator it = this->getNodeFwIterator(); 
       it!=this->getNodesEnd(); ++it) {
    delete *it;
  }
}



/******************************************************************
 * Basic management of nodes
 *****************************************************************/

// Returns a pointer to the first node i.e. the lowest in the tree
Node* Forest::getFirstNode() {
  return(this->nodes_[0]);
}

// Returns a pointer the first node i.e. the lowest in the tree
std::vector<Node*>::iterator Forest::getNodesEnd() {
  return(this->nodes_.end());
}

// Iterator for moving true nodes sorted by height
std::vector<Node*>::iterator Forest::getNodeFwIterator() {
  return(nodes_.begin());
}

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
    this->inc_local_tree_length(node->height_above());
    this->inc_total_tree_length(node->height_above());
  }

  if (lower_child != NULL) {
    node->set_lower_child(lower_child);
    lower_child->set_parent(node);
    this->inc_local_tree_length(lower_child->height_above());
    this->inc_total_tree_length(lower_child->height_above());
  }

  if (higher_child != NULL) {
    node->set_higher_child(higher_child);
    higher_child->set_parent(node);
    this->inc_local_tree_length(higher_child->height_above());
    this->inc_total_tree_length(higher_child->height_above());
  }
}

void Forest::printNodes() {
  for(int i = 0; i < this->countNodes(); ++i) {
    std::cout << "Addr: " << this->nodes()[i] << " | ";
    std::cout << "Height: " << this->nodes()[i]->height() << " | ";
    std::cout << "Parent: " << this->nodes()[i]->parent() << std::endl;
  }
}


void Forest::buildInitialTree() {
  this->createSampleNodes();

  RandomGenerator *rg = this->random_generator();
  std::cout << "Preparing coalescence" << std::endl;
  std::vector<Node*> uncoalesced_nodes = this->nodes();
  double time = 0;
  while (uncoalesced_nodes.size() > 1) {
    int node1, node2;
    int n = uncoalesced_nodes.size();
    double rate = (1.0/(2*this->model().population_size()))*n*(n-1)/2.0;
    time += rg->sampleExpo(rate);
    rg->sampleTwoElements(uncoalesced_nodes.size(), &node1, &node2);
    std::cout << "Coalescing Nodes " << node1 << " and " << node2 << std::endl;

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
  this->set_ultimate_root(uncoalesced_nodes[0]);
  assert(this->checkTree());
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
  assert(local_tree_length_ > 0);
}

void Forest::inc_total_tree_length(const double &by) {
  total_tree_length_ = total_tree_length_ + by;
  assert(total_tree_length_ > 0);
}


TreePoint Forest::samplePoint(bool only_local) {
  //O(#Nodes) implementation
  double length = 0;

  if (only_local) length = this->local_tree_length();
  else            length = this->total_tree_length();

  double point = this->random_generator()->sample() * length;

  Node* above_node;
  double height_above;

  for (std::vector<Node*>::iterator it = nodes_.begin(); it!=nodes_.end(); ++it) {
    if ((*it)->height_above() > point) {
      above_node = *it;
      height_above = point;
      break;
    }
    point -= (*it)->height_above();
  }

  return TreePoint(above_node, height_above);
}


void Forest::sampleNextGenealogy() {
  // Samples a new genealogy, conditional on a recombination occuring

  // Sample the recombination point
  TreePoint rec_point = this->samplePoint();
  std::cout << "Recombination at " << rec_point.height_above() << " above " << rec_point.above_node() << std::endl;

  // Postpone coalesence if not active
  //if (!idx->active()) { ; }


}



/******************************************************************
 * Debugging Utils
 *****************************************************************/

void Forest::createExampleTree() {
  this->nodes_.clear();
  this->set_local_tree_length(0);
  this->set_total_tree_length(0);
  Node* leaf1 = new Node(0);
  Node* leaf2 = new Node(0);
  Node* leaf3 = new Node(0);
  Node* leaf4 = new Node(0);
  this->addNode(leaf1);
  this->addNode(leaf2);
  this->addNode(leaf3);
  this->addNode(leaf4);

  Node* node12 = new Node(1);
  this->addNodeToTree(node12, NULL, leaf1, leaf2);
  
  Node* node34 = new Node(3);
  this->addNodeToTree(node34, NULL, leaf3, leaf4);
  
  Node* root = new Node(10);
  this->addNodeToTree(root, NULL, node12, node34);
  this->set_ultimate_root(root);
}


bool Forest::checkNodesSorted() {
  double cur_height = 0;
  for (std::vector<Node*>::iterator it = this->getNodeFwIterator(); 
       it!=this->getNodesEnd(); ++it) {

    if ((*it)->height() < cur_height) {
      std::cerr << "Nodes not sorted" << std::endl;
      return(0);
    } else {
      cur_height = (*it)->height();
    }

  }
  return(1);
}

bool Forest::checkTree(Node *root) {
  bool sorted = 1;
  if (root == NULL) {
    //Default when called without argument
    root = this->ultimate_root();

    //Also check if nodes are sorted by height in this case
    sorted = this->checkNodesSorted();
  }

  Node* h_child = root->higher_child();
  Node* l_child = root->lower_child();

  if (h_child != NULL && l_child != NULL) {
    if (h_child->height() < l_child->height()) { 
      std::cerr << "Child Nodes in wrong order" << std::endl;
      return 0;
    }
  }
  else if (! (h_child == NULL && l_child == NULL)) { 
    std::cerr << "Node has only one child" << std::endl;
    return 0;
  }
  
  bool child1 = 1;
  if (h_child != NULL) {
    if (h_child->parent() != root) {
      std::cerr << "Node is child of non-parent" << std::endl;
      return 0;
    }
    child1 = checkTree(h_child);
  }

  bool child2 = 1;
  if (l_child != NULL) {
    if (l_child->parent() != root) {
      std::cerr <<"Node is child of non-parent" << std::endl;
      return 0;
    }
    child2 = checkTree(h_child);
  }

  return sorted*child1*child2;
}

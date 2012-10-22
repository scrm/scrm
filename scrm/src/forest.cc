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
  std::vector<Node*>::iterator it;
  for (it = nodes_.begin(); it!=nodes_.end(); ++it) {
    if ((*it)->height() > node->height()) break;
  }
  nodes_.insert(it, node);
  assert(this->checkNodesSorted());
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
  for (int i=0; i < this->sample_size(); i++) {
    Node* node = new Node(0);
    std::cout << "New node : " << node << std::endl;
    this->addNode(node);
    
    if (i == 0) {
      this->set_ultimate_root(node);      
      continue;
    }

    TreePoint new_leaf = TreePoint(node, 0, true);
    this->sampleCoalescences(new_leaf, true);
  }
}
 
void Forest::buildInitialTree_old() {
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

  Node* base_node;
  double height_above;

  for (std::vector<Node*>::iterator it = nodes_.begin(); it!=nodes_.end(); ++it) {
    if ((*it)->height_above() > point) {
      base_node = *it;
      height_above = point;
      break;
    }
    point -= (*it)->height_above();
  }

  return TreePoint(base_node, height_above, true);
}


void Forest::sampleNextGenealogy() {
  // Samples a new genealogy, conditional on a recombination occuring

  // Sample the recombination point
  TreePoint rec_point = TreePoint(new Node(0), 0, true);
  this->addNode(rec_point.base_node());
  std::cout << "Recombination at " << rec_point.relative_height() << " above " << rec_point.base_node() << std::endl;

  // Postpone coalesence if not active
  //if (!idx->active()) { ; }

  this->sampleCoalescences(rec_point);
}

void Forest::sampleCoalescences(TreePoint &start_point, const bool &for_initial_tree) {
  EventIterator events = EventIterator(this, start_point.height());

  //Node* coalescence_root_1 = new Node(start_point.height());
  Node* coalescence_root_1 = start_point.base_node();
  Node* coalescence_root_2 = this->local_root();
  bool  coalescence_1_active = false;
  bool  coalescence_2_active = false;

  double current_time = std::min(start_point.height(), coalescence_root_2->height());

  Event current_event = events.next();
  double std_expo_sample = this->random_generator()->sampleExpo(1);
  double expo_sample = -1;

  while (true) {
    assert(this->checkTree());
    std::cout << "New Event - Time: " << current_time << std::endl;

    coalescence_1_active = current_time >= coalescence_root_1->height();
    coalescence_2_active = current_time >= coalescence_root_2->height();
    if (coalescence_1_active) std::cout << "L1 active" << std::endl;
    if (coalescence_2_active) std::cout << "L2 active" << std::endl;
    assert(coalescence_1_active || coalescence_2_active);

    // If SMC or intiatial tree:
    // Remove branch above recombination point to avoid back coalescence.
    if ( for_initial_tree  || this->model().is_smc_model() ) {
      current_event.removeFromContemporaries(start_point.base_node());
    }

    // Calculate rate of any coalescence occuring
    double rate = calcCoalescenceRate(current_event.contemporaries().size(),
                                      coalescence_1_active + coalescence_2_active);

    double intervall_height = current_event.end_height() - current_event.start_height();
    expo_sample = std_expo_sample / rate;

    if (expo_sample > intervall_height) {
      std::cout << "No coalescence" << std::endl;
      std_expo_sample = (expo_sample - intervall_height) * rate;
    }
    else {
      current_time += expo_sample;
      std::cout << "Coalescence at time " << current_time << std::endl;
      break;
    }

    //We should at least coalescence in the last (almost infinite) intervall
    assert(current_event.end_height() < FLT_MAX);
    
    current_time = current_event.end_height();
    current_event = events.next();
  }
      
  std::cout << "Contemoraries: " << current_event.contemporaries().size() << std::endl;

  TreePoint coal_point;
  if (coalescence_1_active && coalescence_2_active) {
    // DOTO: Modification to include non-local trees
    coal_point = TreePoint(coalescence_root_2, current_time, false);
    coalesNodeIntoTree(coalescence_root_1, coal_point);
  } else {
    coal_point = TreePoint(current_event.getRandomContemporary(), current_time, false);
    //if (coal_point.base_node() != coalescence_root_1) {
    // ARE NODES ALLOWED TO COALESCENCE INTO THEIR OWN BRANCH?
      coalesNodeIntoTree(coalescence_root_1, coal_point);
    //}
  }



}

// Calculates the rate of any coalescence occuring in an intervall with a total of 
// lines_number lines out of which coal_lines_number are not coalescenced.
double Forest::calcCoalescenceRate(int lines_number, int coal_lines_number) {
  assert(lines_number > 0);
  assert(coal_lines_number > 0 && coal_lines_number <= 2);

  double rate = -1;

  if (coal_lines_number == 2)
    rate = ( 2 * lines_number + 1 ) / ( 2.0 * this->model().population_size() );
  else
    rate = lines_number / ( 2.0 * this->model().population_size() );

  assert(rate > 0);
  return(rate);
}

void Forest::coalesNodeIntoTree(Node* coal_node, TreePoint &coal_point) {
  Node *new_parent = new Node(coal_point.height());
  std::cout << "New Node: " << new_parent << std::endl;
  coal_node->set_parent(new_parent);
  
  Node *higher_child, *lower_child;
  if (coal_node->height() >= coal_point.base_node()->height()) {
    higher_child = coal_node;
    lower_child = coal_point.base_node();
  } else {
    higher_child = coal_point.base_node();
    lower_child =  coal_node; 
  }

  //Update above
  new_parent->set_parent(coal_point.base_node()->parent());
  if (new_parent->parent() != NULL) {
    //New_parent is not the new root
    new_parent->parent()->change_child(coal_point.base_node(), new_parent);
  } else {
    this->set_ultimate_root(new_parent);
  }

  //Update childs
  new_parent->set_higher_child(higher_child);
  new_parent->set_lower_child(lower_child);
  higher_child->set_parent(new_parent);
  lower_child->set_parent(new_parent);


  //this->addNode(coal_node);
  this->addNode(new_parent);

  this->printNodes();
  assert(this->checkTree());
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
      std::cerr << root << ": Child Nodes in wrong order" << std::endl;
      std::cerr << root << ": higher child " << h_child << " at " << h_child->height() << std::endl;
      std::cerr << root << ": lower child " << l_child << " at " << l_child->height() << std::endl;
      return 0;
    }
  }
  else if (! (h_child == NULL && l_child == NULL)) { 
    std::cerr << root << ": Has only one child" << std::endl;
    return 0;
  }
  
  bool child1 = 1;
  if (h_child != NULL) {
    if (h_child->parent() != root) {
      std::cerr << h_child << ": is child of non-parent" << std::endl;
      return 0;
    }
    child1 = checkTree(h_child);
  }

  bool child2 = 1;
  if (l_child != NULL) {
    if (l_child->parent() != root) {
      std::cerr << l_child << ": is child of non-parent" << std::endl;
      return 0;
    }
    child2 = checkTree(h_child);
  }

  return sorted*child1*child2;
}

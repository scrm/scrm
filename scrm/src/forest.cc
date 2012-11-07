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
//    delete *it;
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

void Forest::removeNode(Node *node) {
  std::vector<Node*>::iterator it;
  for (it = nodes_.begin(); it!=nodes_.end(); ++it) {
    if ((*it) == node ) break;
  }
  if (it == nodes_.end()) throw std::logic_error("Trying to delete apparently non-existing node");
  nodes_.erase(it);
}

void Forest::moveNode(Node *node, const double &new_height) {
  this->removeNode(node);
  this->addNode(node);
}

void Forest::printNodes() {
  for(int i = 0; i < this->countNodes(); ++i) {
    dout << "Addr: " << this->nodes()[i] << " | ";
    dout << "Height: " << this->nodes()[i]->height() << " | ";
    if (this->nodes()[i]->is_ultimate_root())
      dout << "Parent: " << "NONE" << " | ";
    else
      dout << "Parent: " << this->nodes()[i]->parent() << " | ";
    dout << "h_child: " << this->nodes()[i]->higher_child() << " | ";
    dout << "l_child: " << this->nodes()[i]->lower_child() << std::endl;
  }
}

/******************************************************************
 * Tree Printing
 *****************************************************************/
void Forest::printTree() {
  std::vector<Node*>::reverse_iterator it;
  std::vector<Node*>::iterator cit;
  Node* current_node;
  int lines_left, lines_right, position;
  
  std::vector<Node*> branches;
  for (int i=0; i < this->countNodes(); ++i) branches.push_back(NULL);

  for (it = nodes_.rbegin(); it < nodes_.rend(); ++it) {
    current_node = *it;
    if (current_node->is_fake()) continue;


    lines_left = countLinesLeft(current_node);
    lines_right = countLinesRight(current_node);
    

    if ( current_node->is_root() ) {
      position = countBelowLinesLeft(current_node) + lines_left;
      branches[position] = current_node;
    } else {
       position = 0;
       for (cit = branches.begin(); cit < branches.end(); ++cit) {
         if ( *cit == current_node ) break;
         ++position;
       }
    }

    for (cit = branches.begin(); cit < branches.end(); ++cit) {
      if (*cit == NULL) dout << " ";
      else dout << "|";
    }
    dout << std::endl;
    
    for (int i=0; i < branches.size(); ++i) {
      if ( i < position - lines_left) {
        if ( branches[i] == NULL ) dout << " ";
        else if (areSame(branches[i]->height(), current_node->height())) 
          dout << "+";
        else dout << "|";
      } 
      else if ( i < position) dout << "-";
      else if ( i == position) dout << "+";
      else if ( i <= position + lines_right) dout << "-";
      else {
        if ( branches[i] == NULL ) dout << " ";
        else if (areSame(branches[i]->height(), current_node->height())) 
          dout << "+";
        else dout << "|";
      }
    }
    
    dout << position << " " << lines_left << " " << lines_right;
    ////dout << " " << countLinesLeft(current_node);
    dout << std::endl;

    branches[position - lines_left] = current_node->lower_child();
    branches[position + lines_right] = current_node->higher_child();
    branches[position] = NULL;
    if (current_node ->height() == 0) break;
  }
}

int Forest::countLinesLeft(Node *node) {
  if ( node->lower_child() == NULL ) return 0;
  return ( 1 + countBelowLinesRight(node) );
}

int Forest::countLinesRight(Node *node) {
  if ( node->lower_child() == NULL ) return 0;
  return ( 1 + countBelowLinesLeft(node) );
}

int Forest::countBelowLinesLeft(Node *node) {
  if ( node->lower_child() == NULL ) return 0;
  else return ( countLinesLeft(node->lower_child()) + countBelowLinesLeft(node->lower_child()) );
}

int Forest::countBelowLinesRight(Node *node) {
  if ( node->higher_child() == NULL ) return 0;
  else return ( countLinesRight(node->higher_child()) + countBelowLinesRight(node->higher_child()) );
}


//Cuts all nodes below a point on the tree and moves them into a new tree
void Forest::cut(const TreePoint &cut_point) {
  //The node above the cut_point in the old tree
  Node* parent = cut_point.base_node()->parent();
  assert( parent != NULL );

  //The new end of the old branch after the cut
  Node* new_leaf = new Node(cut_point.height(), false);
  new_leaf->set_parent(parent);
  parent->change_child(cut_point.base_node(), new_leaf);
  this->addNode(new_leaf);

  //The new "root" of the newly formed tree
  Node* new_root = new Node(cut_point.height());
  cut_point.base_node()->set_parent(new_root);
  new_root->set_higher_child(cut_point.base_node());

  this->updateTreeLength();
  assert( this->checkTree() );
}

void Forest::calcTreeLength(double *local_length, double *total_length) {
  *local_length = 0;
  *total_length = 0;

  for (std::vector<Node*>::iterator it = this->getNodeFwIterator(); 
       it!=this->getNodesEnd(); ++it) {
    if ( (*it)->is_fake() || (*it)->is_root() ) continue;
    *total_length += (*it)->height_above();
    if ( (*it)->active() ) *local_length += (*it)->height_above();
  }
}

void Forest::updateTreeLength() {
  double local_length = 0, total_length = 0;
  this->calcTreeLength(&local_length, &total_length);
  this->set_local_tree_length(local_length);
  this->set_total_tree_length(total_length);
}

void Forest::buildInitialTree() {
  dout << "===== BUILDING INTITIAL TREE =====" << std::endl;
  dout << "* creating roots... ";
  this->createRoots();
  dout << "done." << std::endl;

  this->printTree();

  dout << "* adding a node ";
  for (int i=1; i < this->sample_size(); i++) {
    //Create a new sepearte little tree of and at height zero
    Node* new_leaf = new Node(0);
    Node* new_root = new Node(0);
    dout << "(" << new_leaf << ")" << std::endl;
    this->addNode(new_leaf);
    this->addNode(new_root);
    new_leaf->set_parent(new_root);
    new_root->set_higher_child(new_leaf);
    dout << "* * staring coalesces" << std::endl;
    
    //Coales the seperate tree into the main tree
    this->sampleCoalescences(new_root, true);
  }
}


void Forest::buildInitialTree_old() {
  RandomGenerator *rg = this->random_generator();
  dout << "Preparing coalescence" << std::endl;
  std::vector<Node*> uncoalesced_nodes = this->nodes();
  double time = 0;
  while (uncoalesced_nodes.size() > 1) {
    int node1, node2;
    int n = uncoalesced_nodes.size();
    double rate = (1.0/(2*this->model().population_size()))*n*(n-1)/2.0;
    time += rg->sampleExpo(rate);
    rg->sampleTwoElements(uncoalesced_nodes.size(), &node1, &node2);
    dout << "Coalescing Nodes " << node1 << " and " << node2 << std::endl;

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

void Forest::inc_tree_length(const double &by, const bool &active) {
  total_tree_length_ += by;
  if (active) local_tree_length_ += by;
  assert(local_tree_length_ > 0);
  assert(total_tree_length_ > 0);
}

TreePoint Forest::samplePoint(bool only_local) {
  //O(#Nodes) implementation
  double length = 0;

  if (only_local) length = this->local_tree_length();
  else            length = this->total_tree_length();

  double point = this->random_generator()->sample() * length;

  Node* base_node = NULL;
  double height_above = -1;

  for (std::vector<Node*>::iterator it = nodes_.begin(); it!=nodes_.end(); ++it) {
    if ((*it)->height_above() > point) {
      base_node = *it;
      height_above = point;
      break;
    }
    point -= (*it)->height_above();
  }

  assert(base_node != NULL);
  assert(height_above >= 0);
  return TreePoint(base_node, height_above, true);
}


void Forest::sampleNextGenealogy() {
  // Samples a new genealogy, conditional on a recombination occuring

  // Sample the recombination point
  TreePoint rec_point = this->samplePoint();

  dout << "Recombination at " << rec_point.relative_height() << " ";
  dout << "above " << rec_point.base_node() << std::endl;

  this->cut(rec_point);

  // Postpone coalesence if not active
  //if (!idx->active()) { ; }

  //this->sampleCoalescences(rec_point);
}

void Forest::sampleCoalescences(Node *start_node, const bool &for_initial_tree) {
  EventIterator events = EventIterator(this, start_node->height());

  assert(start_node->is_root());

  //Node* coalescence_root_1 = new Node(start_point.height());
  Node* coalescence_root_1 = start_node;
  Node* coalescence_root_2 = this->local_root();
  bool  coalescence_1_active = false;
  bool  coalescence_2_active = false;

  double current_time = std::min(start_node->height(), coalescence_root_2->height());

  Event current_event = events.next();
  double std_expo_sample = this->random_generator()->sampleExpo(1);
  double expo_sample = -1;

  while (true) {
    assert(this->checkTree());
    dout << "* * new event (time " << current_time << ")" << std::endl;

    coalescence_1_active = current_time >= coalescence_root_1->height();
    coalescence_2_active = current_time >= coalescence_root_2->height();
    assert(coalescence_1_active || coalescence_2_active);

    // If SMC or intiatial tree:
    // Remove branch above recombination point to avoid back coalescence.
    if ( for_initial_tree  || this->model().is_smc_model() ) {
      current_event.removeFromContemporaries(start_node);
    }

    // Calculate rate of any coalescence occuring
    double rate = calcCoalescenceRate(current_event.contemporaries().size(),
                                      coalescence_1_active + coalescence_2_active);
    dout << "* * * rate: " << rate << std::endl;

    double intervall_height = current_event.end_height() - current_event.start_height();
    expo_sample = std_expo_sample / rate;

    if (expo_sample > intervall_height) {
      std_expo_sample = (expo_sample - intervall_height) * rate;
    }
    else {
      current_time += expo_sample;
      dout << "* * * coalescence" << std::endl;
      break;
    }

    //We should at least coalescence in the last (almost infinite) intervall
    assert(current_event.end_height() < FLT_MAX);
    
    current_time = current_event.end_height();
    current_event = events.next();
  }
      
  dout << "* * coalescence at time " << current_time << std::endl;
  dout << "* * #contemoraries: " << current_event.contemporaries().size() << std::endl;

  TreePoint coal_point;
  if (coalescence_1_active && coalescence_2_active) {
    // DOTO: Modification to include non-local trees
    coal_point = TreePoint(coalescence_root_2, current_time, false);
    dout << "* * into: " << coal_point.base_node() << std::endl;
    coalesNodeIntoTree(coalescence_root_1, coal_point);
  } else {
    coal_point = TreePoint(current_event.getRandomContemporary(), current_time, false);
    dout << "* * into: " << coal_point.base_node() << std::endl;
    coalesNodeIntoTree(coalescence_root_1, coal_point);
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

void Forest::coalesNodeIntoTree(Node* coal_node, const TreePoint &coal_point) {
  //Move coal_node up to its new position in the tree
  coal_node->set_height(coal_point.height());
  this->moveNode(coal_node, coal_point.height());

  //Update the coal_node
  coal_node->change_child(NULL, coal_point.base_node());
  coal_node->set_parent(coal_point.base_node()->parent());

  //Update the new child 
  coal_point.base_node()->set_parent(coal_node);

  //Update the parent
  coal_node->parent()->change_child(coal_point.base_node(), coal_node);

  //Optimize: Live tracking of tree length?
  double local_length = 0;
  double total_length = 0;
  for (std::vector<Node*>::iterator it = this->getNodeFwIterator(); 
       it!=this->getNodesEnd(); ++it) {
    if ( (*it)->is_fake() || (*it)->is_root() ) continue;
    total_length += (*it)->height_above();
    if ( (*it)->active() ) local_length += (*it)->height_above();
  }
  this->set_local_tree_length(local_length);
  this->set_total_tree_length(total_length);

  this->printTree();
  assert(this->checkTree());
}


/******************************************************************
 * Manage roots and fake tree
 *****************************************************************/

// Creates the very first node of the tree and the ultimate root above
void Forest::createRoots() {
  //Create the first node of the tree (also is the root at this point)
  Node* first_node = new Node(0);
  this->addNode(first_node);

  //Create the ultimate root
  Node* ultimate_root = new Node(FLT_MAX);
  ultimate_root->deactivate();
  this->addNode(ultimate_root);
  this->set_ultimate_root(ultimate_root);

  first_node->set_parent(ultimate_root);
  ultimate_root->set_lower_child(first_node);
  
  assert(this->checkTree());
}

void Forest::registerNonLocalRoot(Node* node) {};
void Forest::unregisterNonLocalRoot(Node* node) {};


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
  
  //Build the fake tree above
  Node* ultimate_root = new Node(FLT_MAX, false);
  this->set_ultimate_root(ultimate_root);
  this->addNode(ultimate_root);

  ultimate_root->set_lower_child(root);
  root->set_parent(ultimate_root);

  assert( this->checkTree() );
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
    this->inc_tree_length(node->height_above(), node->active());
  }

  if (lower_child != NULL) {
    node->set_lower_child(lower_child);
    lower_child->set_parent(node);
    this->inc_tree_length(lower_child->height_above(), lower_child->active());
  }

  if (higher_child != NULL) {
    node->set_higher_child(higher_child);
    higher_child->set_parent(node);
    this->inc_tree_length(higher_child->height_above(), higher_child->active());
  }
}

bool Forest::checkNodesSorted() {
  double cur_height = 0;
  for (std::vector<Node*>::iterator it = this->getNodeFwIterator(); 
       it!=this->getNodesEnd(); ++it) {

    if ((*it)->height() < cur_height) {
      dout << "Error: Nodes not sorted" << std::endl;
      return(0);
    } else {
      cur_height = (*it)->height();
    }

  }
  return(1);
}

bool Forest::checkTreeLength() {
  double local_length = 0;
  double total_length = 0;
  this->calcTreeLength(&local_length, &total_length);

  if ( !areSame(total_length, total_tree_length()) ) {
    dout << "Error: total tree length is " << this->total_tree_length() << " ";
    dout << "but should be " << total_length << std::endl;
    return(0);
  }
  if ( !areSame(local_length, local_tree_length()) ) {
    dout << "Error: local tree length is " << this->local_tree_length() << " ";
    dout << "but should be " << local_length << std::endl;
    return(0);
  }
  return(1);
}

bool Forest::checkTree(Node *root) {
  if (root == NULL) {
    //Default when called without argument
    root = this->ultimate_root();

    //Also check if nodes are sorted by height in this case
    assert(this->checkNodesSorted());
    assert(this->checkTreeLength());
  }
  assert( root != NULL);

  Node* h_child = root->higher_child();
  Node* l_child = root->lower_child();

  //Do we have two childs?
  if (h_child != NULL && l_child != NULL) {
    if (h_child->height() < l_child->height()) { 
      dout << root << ": Child Nodes in wrong order" << std::endl;
      dout << root << ": higher child " << h_child << " at " << h_child->height() << std::endl;
      dout << root << ": lower child " << l_child << " at " << l_child->height() << std::endl;
      return 0;
    }
  }
  //Do we have only one child?
  else if ( !(h_child == NULL && l_child == NULL) ) {
    //This is only allowed for root nodes
    if ( !(root->is_root() || root->is_ultimate_root()) ) { 
      dout << root << ": Has only one child" << std::endl;
      return 0;
    }
  }
  
  bool child1 = 1;
  if (h_child != NULL) {
    if (h_child->parent() != root) {
      dout << h_child << ": is child of non-parent" << std::endl;
      return 0;
    }
    child1 = checkTree(h_child);
  }

  bool child2 = 1;
  if (l_child != NULL) {
    if (l_child->parent() != root) {
      dout << l_child << ": is child of non-parent" << std::endl;
      return 0;
    }
    child2 = checkTree(l_child);
  }

  return child1*child2;
}

bool areSame(double a, double b) {
      return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
      //return std::fabs(a - b) < 0.0001;
}

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
    if ((*it)->height() >= node->height()) break;
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
  dout << "* * New leaf of local tree: " << new_leaf << std::endl;

  //The new "root" of the newly formed tree
  Node* new_root = new Node(cut_point.height(), 
                            cut_point.base_node()->active());
  cut_point.base_node()->set_parent(new_root);
  new_root->set_lower_child(cut_point.base_node());
  //Inefficient to to this for active nodes
  this->addNode(new_root);
  this->registerNonLocalRoot(new_root);
  assert(this->printNodes());
  assert(this->printTree());
  dout << "* * New root of subtree: " << new_root << std::endl;

  this->updateTreeLength();
  assert( this->checkTree() );
  dout << "* * Done" << std::endl;
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

  assert(printTree());

  for (int i=1; i < this->sample_size(); i++) {
    dout << "* adding node ";
    //Create a new sepearte little tree of and at height zero
    Node* new_leaf = new Node(0);
    Node* new_root = new Node(0);
    dout << "(" << new_leaf << ")" << std::endl;
    this->addNode(new_leaf);
    this->addNode(new_root);
    registerNonLocalRoot(new_root);
    new_leaf->set_parent(new_root);
    new_root->set_higher_child(new_leaf);
    dout << "* * staring coalesces" << std::endl;
    
    //Coales the seperate tree into the main tree
    this->sampleCoalescences(new_root, false);
    dout << "* * Tree:" << std::endl;
    assert(this->printTree());
  }
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
    if ((*it)->is_root()) continue;
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
  dout << "===== BUILDING NEXT GENEALOGY =====" << std::endl;
  // Samples a new genealogy, conditional on a recombination occurring

  // Sample the recombination point
  TreePoint rec_point = this->samplePoint();

  dout << "* Recombination at " << rec_point.relative_height() << " ";
  dout << "above " << rec_point.base_node() << std::endl;

  dout << "* Cutting subtree below recombination " << std::endl;
  this->cut(rec_point);

  //assert(this->printTree());

  // Postpone coalescence if not active
  if (!rec_point.base_node()->active()) {
    dout << "* Not on local tree; Postponing coalescence" << std::endl; 
    dout << "* Tree:" << std::endl; 
    assert(this->printTree());
    return;
  }

  dout << "* Starting coalescence" << std::endl;
  this->sampleCoalescences(rec_point.base_node()->parent());
  //assert(this->printTree());
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
  double std_expo_sample = -1;
  double expo_sample = -1;

  while (true) {
    Event current_event = events.next();
    double std_expo_sample = this->random_generator()->sampleExpo(1);

    while (true) {
      assert(this->checkTree());
      dout << "* * new event (time " << current_time << ")" << std::endl;

      coalescence_1_active = current_time >= coalescence_root_1->height();
      coalescence_2_active = current_time >= coalescence_root_2->height();
      assert(coalescence_1_active || coalescence_2_active);
      dout << "* * Active: ";
      if (coalescence_1_active) dout << coalescence_root_1 << " ";
      if (coalescence_2_active) dout << coalescence_root_2 << " ";
      dout << std::endl;

      // If SMC or intiatial tree:
      // Remove branch above recombination point to avoid back coalescence.
      if ( for_initial_tree  || this->model().is_smc_model() ) {
        dout << "* * removing branch above coalescing node" << std::endl;
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
 
    Node *coalescing_root = NULL, *coalescence_target = NULL;
    
    if (coalescence_1_active && coalescence_2_active) {
      // Select which root coaleses
      if (current_event.contemporaries().size() == 0 ||
          this->random_generator()->sampleInt(2) == 0) {
        coalescing_root = coalescence_root_1;
        // root_1 can also coales into root_2, while the opposit is not possible.
        // (Or the pair would be overrepresented)
        int target = this->random_generator()->sampleInt(current_event.contemporaries().size() + 1);
        if (target == 0) coalescence_target = coalescence_root_2;
        else coalescence_target = current_event.contemporaries()[target-1];
      } else {
        coalescing_root = coalescence_root_2;
        coalescence_target = current_event.getRandomContemporary();
      }
    } else {
      // One active coalescence
      assert(coalescence_1_active || coalescence_2_active);
      if (coalescence_1_active) coalescing_root = coalescence_root_1;
      else                      coalescing_root = coalescence_root_2;
      coalescence_target = current_event.getRandomContemporary();
    }

    assert(coalescing_root != NULL);
    assert(coalescence_target != NULL);
    assert(coalescing_root != coalescence_target);
    dout << "* * " << coalescing_root << " >> " << coalescence_target << std::endl;

    // Implement the coalescence
    TreePoint coal_point = TreePoint(coalescence_target, current_time, false);
    dout << "* * into: " << coal_point.base_node() << std::endl;
    dout << "* * Updating tree" << std::endl;
    coalesNodeIntoTree(coalescing_root, coal_point);

    // root 1 and 2 coalesed together? => we are done
    if (coalescing_root == coalescence_root_1 && coalescence_target == coalescence_root_2) {
      dout << "* * pairwise coalescence" << std::endl;
      break;
    }

    // root 2 coalesed? => move upwards until we hid another root, and mark as
    // active
    if (coalescing_root == coalescence_root_2) {
      while (!coalescence_root_2->is_root()) {
        coalescence_root_2 = coalescence_root_2->parent();
        coalescence_root_2->activate();
        dout << "* * Activating node " << coalescence_root_2 << std::endl;
      }
      dout << "* * New second root: " << coalescence_root_2 << std::endl;
    }

    // root 1 coalesed? follow branch up until
    // - we hit ancestral branch => Done
    // - we hit a root => continue coalescence
    if (coalescing_root == coalescence_root_1) {
      while ( !(coalescence_root_1->is_root()) ) {
        coalescence_root_1 = coalescence_root_1->parent();
        if (coalescence_root_1->active()) break;
        dout << "* * Activating node " << coalescence_root_1 << std::endl;
        coalescence_root_1->activate();
      }
      if (coalescence_root_1->active()) {
        dout << "* * Hit active node " << coalescence_root_1 << std::endl;
        break;
      }
      assert(coalescence_root_1->is_root());
    }

    if (coalescence_root_1 == coalescence_root_2) {
      dout << "* We coalesced, apparently." << std::endl;
      break;
    }

    current_time = std::max(current_time, std::min(coalescence_root_1->height(), 
                                                   coalescence_root_2->height()));

  }

  dout << "* * Tree:" << std::endl;
  assert(this->printTree());
}

// Calculates the rate of any coalescence occuring in an intervall with a total of 
// lines_number lines out of which coal_lines_number are not coalescenced.
double Forest::calcCoalescenceRate(int lines_number, int coal_lines_number) {
  assert(lines_number >= 0);
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
  dout << "* * * Moving root of subtree to it new position" << std::endl;
  coal_node->set_height(coal_point.height());

  this->moveNode(coal_node, coal_point.height());

  //Update the coal_node
  dout << "* * * Updating its structure" << std::endl;
  coal_node->change_child(NULL, coal_point.base_node());
  coal_node->set_parent(coal_point.base_node()->parent());
  unregisterNonLocalRoot(coal_node);

  //Update the new child 
  coal_point.base_node()->set_parent(coal_node);

  //Update the parent
  coal_node->parent()->change_child(coal_point.base_node(), coal_node);

  //Optimize: Live tracking of tree length?
  this->updateTreeLength();

  assert(this->checkTree());
  dout << "* * * Done" << std::endl;
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

void Forest::registerNonLocalRoot(Node* node, Node *position) {
  if (position == NULL) position = this->ultimate_root();
  //Do we have a free place to add the node?
  if (position->lower_child() == NULL) {
    position->set_lower_child(node);
    node->set_parent(position);
    return;
  }

  //If not, we need to have another fake node immediately below us 
  Node* next_fake_node = position->higher_child();
  if (next_fake_node == NULL) {
    next_fake_node = new Node(FLT_MAX);
    next_fake_node->set_parent(position);
    position->set_higher_child(next_fake_node);
    this->addNode(next_fake_node);
  }

  registerNonLocalRoot(node, next_fake_node);
};

void Forest::unregisterNonLocalRoot(Node* node, Node* position) {
  if (position == NULL) position = this->ultimate_root();
  if (position->lower_child() == node) position->set_lower_child(NULL);
  else {
    if (position->higher_child() == NULL)
      throw std::logic_error("Cannont unregister root: Not found");
    unregisterNonLocalRoot(node, position->higher_child());
  } 
};




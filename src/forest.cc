#include "forest.h"


/******************************************************************
 * Constructors & Initialization
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

  this->nodes_ = new NodeContainer();
  this->set_model(model);
  this->set_random_generator(rg);
  this->set_ultimate_root(ultimate_root);
  this->set_local_tree_length(local_tree_length);
  this->set_total_tree_length(total_tree_length);
  this->set_current_base(1);
}



/******************************************************************
 * Destructor
 *****************************************************************/

Forest::~Forest() { 
  delete nodes_;
}

//Cuts all nodes below a point on the tree and moves them into a new tree
void Forest::cut(const TreePoint &cut_point) {
  //The node above the cut_point in the old tree
  Node* parent = cut_point.base_node()->parent();
  assert( parent != NULL );

  //The new end of the old branch after the cut
  Node* new_leaf = new Node(cut_point.height(), false, current_base(), 0, 0);
  new_leaf->set_parent(parent);
  parent->change_child(cut_point.base_node(), new_leaf);
  nodes()->add(new_leaf);
  updateAbove(parent);
  dout << "* * New leaf of local tree: " << new_leaf << std::endl;
  
  //The new "root" of the newly formed tree
  Node* new_root = new Node(cut_point.height(), 
                            cut_point.base_node()->active());
  cut_point.base_node()->set_parent(new_root);
  new_root->set_lower_child(cut_point.base_node());
  //Inefficient to to this for active nodes
  nodes()->add(new_root);
  this->registerNonLocalRoot(new_root);
  assert(this->printTree());
  dout << "* * New root of subtree: " << new_root << std::endl;

  assert( this->checkTree() );
  dout << "* * Done" << std::endl;
}

void Forest::updateAbove(Node* node, bool above_local_root, bool recursive) {
  dout << "* * * Updating: " << node << " fastforward: " << above_local_root << std::endl;
  // Fast forward above local root because this part is non-local
  if (above_local_root) {
    node->deactivate(current_base());
    node->set_samples_below(this->sample_size());
    if ( node->is_root() ) {
      set_primary_root(node);
      return;
    }
    if ( recursive ) updateAbove(node->parent(), true);
  }

  // Calculate new values for samples_below and length_below for the current
  // node
  Node *l_child = node->lower_child();
  Node *h_child = node->higher_child();

  size_t samples_below = 1;
  if (l_child != NULL) samples_below = l_child->samples_below();
  if (h_child != NULL) samples_below += h_child->samples_below();

  double length_below = 0;
  if (l_child != NULL) {
    length_below += l_child->length_below();
    if (l_child->active()) length_below += l_child->height_above();
  }

  if (h_child != NULL) {
    length_below += h_child->length_below();
    if (h_child->active()) length_below += h_child->height_above();
  }
  
  // If nothing changed, we also don't need to update the tree further above.
  if (recursive &&
      samples_below == node->samples_below() && 
      areSame(length_below, node->length_below()) ) return;
  
  // Update the node
  if ( samples_below == 0 ) node->deactivate(current_base());
  else if ( samples_below < this->sample_size() ) node->activate();
  node->set_samples_below(samples_below);
  node->set_length_below(length_below);

  // Check if we are above or equal the local root
  if ( samples_below == this->sample_size() ) {
    node->deactivate(current_base());
    // Are we the local root?
    if (l_child->samples_below() > 0 && h_child->samples_below() > 0) {
      set_local_root(node);
    }
    if ( node->is_root() ) set_primary_root(node);
    above_local_root = true;
  }

  // Go further up if possible
  if ( recursive && !node->is_root() ) {
    updateAbove(node->parent(), above_local_root);
  }
}


// Recursively finds the largest subtree containing "node" which only
// has deactivated leafs. Makes sure that all nodes on this tree
// are deactivated.
/* void Forest::deactivateSubtree(Node* node) {
  if ( node->is_root() ) return;
  if ( !node->higher_child()->active() && !node->lower_child()->active() ) {
    node->deactivate();
    deactivateSubtree(node->parent());
  }
}*/

void Forest::calcTreeLength(double *local_length, double *total_length) const {
  *local_length = 0;
  *total_length = 0;

  for (ConstNodeIterator it = nodes()->iterator(); it.good(); ++it) {
    if ( (*it)->is_fake() || (*it)->is_root() ) continue;
    *total_length += (*it)->height_above();
    if ( (*it)->active() ) *local_length += (*it)->height_above();
  }
}

void Forest::updateTreeLength() {
  assert( false );
  double local_length = 0, total_length = 0;
  this->calcTreeLength(&local_length, &total_length);
//  this->set_local_tree_length(local_length);
//  this->set_total_tree_length(total_length);
}

void Forest::buildInitialTree() {
  dout << "===== BUILDING INITIAL TREE =====" << std::endl;
  dout << "* creating roots... ";
  this->createRoots();
  dout << "done." << std::endl;
  
  for (int i=1; i < this->model().sample_size(); i++) {
    this->set_sample_size(i+1);
    dout << "* adding node ";
    //Create a new separate little tree of and at height zero
    Node* new_leaf = new Node(0, true, 0, 1);
    Node* new_root = new Node(0);
    dout << "(" << new_leaf << ")" << std::endl;
    nodes()->add(new_leaf);
    nodes()->add(new_root);
    registerNonLocalRoot(new_root);
    new_leaf->set_parent(new_root);
    new_root->set_higher_child(new_leaf);
    dout << "* * staring coalesces" << std::endl;
    
    //Coalesces the separate tree into the main tree
    this->sampleCoalescences(new_root, false);
    dout << "* * Tree:" << std::endl;
    assert(this->printNodes());
    assert(this->printTree());
  }
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

  for (NodeIterator it = nodes()->iterator(); it.good(); ++it) {
    if ( (*it)->is_root() ) continue;
    if ( only_local && !(*it)->active() ) continue;

    if ( (*it)->height_above() > point ) {
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
  dout << "Sequence position: " << this->current_base() << std::endl;
  // Samples a new genealogy, conditional on a recombination occurring

  // Sample the recombination point
  TreePoint rec_point = this->samplePoint(true);

  dout << "* Recombination at height " << rec_point.height() << " ";
  dout << "(above " << rec_point.base_node() << ")"<< std::endl;

  dout << "* Cutting subtree below recombination " << std::endl;
  this->cut(rec_point);

  //assert(this->printTree());

  // Postpone coalescence if not active
  if (!rec_point.base_node()->active()) {
    dout << "* Not on local tree; Postponing coalescence" << std::endl; 
    return;
  }

  dout << "* Starting coalescence" << std::endl;
  this->sampleCoalescences(rec_point.base_node()->parent());
  assert(this->printNodes());
  assert(this->printTree());
  assert(this->checkLeafsOnLocalTree());
}

void Forest::sampleCoalescences(Node *start_node, const bool &for_initial_tree) {
  assert(start_node->is_root());

  //Node* coalescence_root_1 = new Node(start_point.height());
  Node* coalescence_root_1 = start_node;
  Node* coalescence_root_2 = this->primary_root();
  bool  coalescence_1_active = false;
  bool  coalescence_2_active = false;

  double current_time = std::min(start_node->height(), coalescence_root_2->height());
  double std_expo_sample = -1;
  double expo_sample = -1;

  while (true) {
    EventIterator events = EventIterator(this, current_time);
    Event current_event = events.next();
    std_expo_sample = this->random_generator()->sampleExpo(1);

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

      // If SMC or initial tree:
      // Remove branch above recombination point to avoid back coalescence.
      if ( for_initial_tree  || this->model().is_smc_model() ) {
        dout << "* * removing branch above coalescing node" << std::endl;
        current_event.removeFromContemporaries(start_node);
      }

      // Calculate rate of any coalescence occurring
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

      //We should at least coalescence in the last (almost infinite) interval
      assert(current_event.end_height() < FLT_MAX);
      
      current_time = current_event.end_height();
      current_event = events.next();
    }
        
    dout << "* * coalescence at time " << current_time << std::endl;
    dout << "* * #contemporaries: " << current_event.contemporaries().size() << std::endl;
 
    Node *coalescing_root = NULL, *coalescence_target = NULL;
    
    if (coalescence_1_active && coalescence_2_active) {
      // Select which root coalesces
      if (current_event.contemporaries().size() == 0 ||
          this->random_generator()->sampleInt(2) == 0) {
        coalescing_root = coalescence_root_1;
        // root_1 can also coalesces into root_2, while the opposite is not possible.
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
      dout << "One active" << std::endl;
      assert(coalescence_1_active || coalescence_2_active);
      if (coalescence_1_active) coalescing_root = coalescence_root_1;
      else                      coalescing_root = coalescence_root_2;
      coalescence_target = current_event.getRandomContemporary();
      dout << "Target: " << coalescence_target << std::endl;
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

    // Root 1 and 2 coalesced together? => we are done
    if (coalescing_root == coalescence_root_1 && coalescence_target == coalescence_root_2) {
      dout << "* * Both roots coalesced together" << std::endl;
      updateAbove(local_root(), false, false);    // To active the old local root
      updateAbove(coalescing_root, false, false); // Update the new node
      break;
    }

    // Root 2 coalesced? Move upwards until we hid another root
    if (coalescing_root == coalescence_root_2) {
      dout << "* * Root 2 coalesced" << std::endl;
      coalescence_root_2 = moveUpwardsInTree(coalescence_root_2);
      coalescence_root_2->activate();
      dout << "* * New Root 2: " << coalescence_root_2 << std::endl;
    }

    // Root 1 coalesced? Follow branch up until
    // - we hit a local branch    => Done
    // - we hit the primary root  => Done 
    // - we hit another root      => continue coalescence
    if (coalescing_root == coalescence_root_1) {
      dout << "* * Root 1 coalesced" << std::endl;
      coalescence_root_1 = moveUpwardsInTree(coalescence_root_1);

      // If it is active, then we hit an active node
      if ( coalescence_root_1->active() || coalescence_root_1 == primary_root() ) {
        dout << "* * * is active or primary root" << std::endl;
        break;
      }

      // Else also active the current node and continue coalescing
      dout << "* * * is another root" << std::endl;
      coalescence_root_1->activate();
      dout << "* * New Root 1: " << coalescence_root_1 << std::endl;
    }

    // Maybe both hit the same third tree...
    if (coalescence_root_1 == coalescence_root_2) {
      dout << "* We coalesced, apparently." << std::endl;
      break;
    }

    current_time = std::max(current_time, std::min(coalescence_root_1->height(), 
                                                   coalescence_root_2->height()));

  }
  dout << "* * Coalescence done" << std::endl;
}

// Starting from "node", this moves upwards in the tree until it hits a
// root or an local node. Updates the nodes on its way up.
// The final root or local node is returned.
Node* Forest::moveUpwardsInTree(Node* node) {
  while ( !node->is_root() ) {
    updateAbove(node, false, false); // Updates only "node"
    node = node->parent();
    dout << "* * * Moving upwards: " << node << std::endl;        
    if (node->active()) {
      updateAbove(node, false); // Updates the rest of the tree
      return node;
    }
  }

  updateAbove(node, false, false); // Updates the found root (only)
  return node;
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

  this->nodes()->move(coal_node, coal_point.height());

  //Update the coal_node
  dout << "* * * Updating its structure" << std::endl;
  coal_node->change_child(NULL, coal_point.base_node());
  coal_node->set_parent(coal_point.base_node()->parent());
  unregisterNonLocalRoot(coal_node);

  //Update the new child 
  coal_point.base_node()->set_parent(coal_node);

  //Update the parent
  coal_node->parent()->change_child(coal_point.base_node(), coal_node);

  //updateAbove(coal_node);

  assert(this->checkTree());
  dout << "* * * Done" << std::endl;
}


/******************************************************************
 * Manage roots and fake tree
 *****************************************************************/

// Creates the very first node of the tree and the ultimate root above
void Forest::createRoots() {
  //Create the first node of the tree (also is the root at this point)
  Node* first_node = new Node(0, true, 0, 1);
  this->nodes()->add(first_node);
  this->set_local_root(first_node);
  this->set_primary_root(first_node);

  //Create the ultimate root
  Node* ultimate_root = new Node(FLT_MAX);
  ultimate_root->deactivate(current_base());
  this->nodes()->add(ultimate_root);
  this->set_ultimate_root(ultimate_root);

  first_node->set_parent(ultimate_root);
  ultimate_root->set_lower_child(first_node);
 
  assert(this->nodes()->size() == 2); 
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
    this->nodes()->add(next_fake_node);
  }

  registerNonLocalRoot(node, next_fake_node);
};

void Forest::unregisterNonLocalRoot(Node* node, Node* position) {
  if (position == NULL) position = this->ultimate_root();
  if (position->lower_child() == node) position->set_lower_child(NULL);
  else if (position->higher_child() == node) position->set_higher_child(NULL);
  else {
    if (position->higher_child() == NULL)
      throw std::logic_error("Cannont unregister root: Not found");
    unregisterNonLocalRoot(node, position->higher_child());
  } 
};




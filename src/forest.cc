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
                        Node* ultimate_root) {

  this->nodes_ = new NodeContainer();
  this->set_model(model);
  this->set_random_generator(rg);
  this->set_ultimate_root(ultimate_root);
  this->set_current_base(1);
  this->expo_sample_ = -1;
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
  // Update all local nodes above.
  updateAbove(parent, false, true, true);
  //this->postpone_update(parent);
  dout << "* * New leaf of local tree: " << new_leaf << std::endl;
  
  //The new "root" of the newly formed tree
  Node* new_root = new Node(cut_point.height(), 
                            cut_point.base_node()->local());
  cut_point.base_node()->set_parent(new_root);
  new_root->set_lower_child(cut_point.base_node());
  //Inefficient to to this for local nodes
  nodes()->add(new_root);
  this->registerNonLocalRoot(new_root);
  assert(this->printTree());
  dout << "* * New root of subtree: " << new_root << std::endl;

  assert( this->checkTree() );
  dout << "* * Done" << std::endl;
}

void Forest::updateAbove(Node* node, bool above_local_root, bool recursive, bool local_only) {
  if (local_only && !node->local()) return;
  dout << "* * * Updating: " << node << " local: " << node->local()
       << " fastforward: " << above_local_root << std::endl;

  // Fast forward above local root because this part is non-local
  if (above_local_root) {
    if (node->local()) node->make_nonlocal(current_base());
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

  size_t samples_below = node->in_sample();
  if (l_child != NULL) samples_below = l_child->samples_below();
  if (h_child != NULL) samples_below += h_child->samples_below();

  double length_below = 0;
  if (l_child != NULL) {
    length_below += l_child->length_below();
    if (l_child->local()) length_below += l_child->height_above();
  }

  if (h_child != NULL) {
    length_below += h_child->length_below();
    if (h_child->local()) length_below += h_child->height_above();
  }
  
  // If nothing changed, we also don't need to update the tree further above.
  if (recursive &&
      samples_below == node->samples_below() && 
      areSame(length_below, node->length_below()) ) return;
  
  // Update the node
  if ( samples_below == 0 ) node->make_nonlocal(current_base());
  else if ( samples_below < this->sample_size() ) node->make_local();
  node->set_samples_below(samples_below);
  node->set_length_below(length_below);

  // Check if we are above or equal the local root
  if ( samples_below == this->sample_size() ) {
    node->make_nonlocal(current_base());
    // Are we the local root?
    if (l_child->samples_below() > 0 && h_child->samples_below() > 0) {
      dout << "* * * is local_root" << std::endl;
      set_local_root(node);
    }
    if ( node->is_root() ) set_primary_root(node);
    above_local_root = true;
  }
  
  // Go further up if possible
  if ( recursive && !node->is_root() ) {
    updateAbove(node->parent(), above_local_root, recursive, local_only);
  }
}


void Forest::buildInitialTree() {
  dout << "===== BUILDING INITIAL TREE =====" << std::endl;
  dout << "* Adding first node... ";
  Node* first_node = new Node(0, true, 0, 1);
  this->nodes()->add(first_node);
  this->set_local_root(first_node);
  this->set_primary_root(first_node);
  dout << "done." << std::endl;
  
  for (int i=1; i < this->model().sample_size(); i++) {
    this->set_sample_size(i+1);
    dout << "* adding node ";
    //Create a new separate little tree of and at height zero
    Node* new_leaf = new Node(0, true, 0, 1);
    dout << "(" << new_leaf << ")" << std::endl;
    nodes()->add(new_leaf);
    dout << "* staring coalesces" << std::endl;
    
    //Coalesces the separate tree into the main tree
    this->sampleCoalescences(new_leaf, false);
    dout << "* * Tree:" << std::endl;
    assert(this->printNodes());
    assert(this->printTree());
  }
}

// uniformly samples a Point on the local tree
// iterative log(#nodes) implementation
// Distribution checked.
TreePoint Forest::samplePoint(Node* node, double length_left) {
  if (node == NULL) {
    // Called without arguments => initialization
    node = this->local_root();
    length_left = this->random_generator()->sample() * node->length_below();
  }

  assert( length_left >= 0 );
  assert( length_left < (node->length_below() + node->height_above()) );
  
  if ( length_left < node->height_above() ) {
    return TreePoint(node, length_left, true);
  }

  length_left = length_left - node->height_above();
  assert( length_left >= 0 );
  double tmp = node->lower_child()->height_above() + node->lower_child()->length_below();
  if ( length_left <= tmp )
    return samplePoint(node->lower_child(), length_left);
  else 
    return samplePoint(node->higher_child(), length_left - tmp);
}


void Forest::sampleNextGenealogy() {
  dout << "===== BUILDING NEXT GENEALOGY =====" << std::endl;
  dout << "Sequence position: " << this->current_base() << std::endl;
  // Samples a new genealogy, conditional on a recombination occurring

  // Sample the recombination point
  TreePoint rec_point = this->samplePoint();

  dout << "* Recombination at height " << rec_point.height() << " ";
  dout << "(above " << rec_point.base_node() << ")"<< std::endl;

  dout << "* Cutting subtree below recombination " << std::endl;
  this->cut(rec_point);

  // Postpone coalescence if not local
  if (!rec_point.base_node()->local()) {
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
  // We can have one or active local nodes: If the coalescing node passes the
  // local root, it also starts a coalescence.
  Node* active_node_1 = start_node;
  Node* active_node_2 = this->local_root();

  // States: Each (branch above a) node can either be in state
  // - 0 = off (the other coalescence has not reached it yet) or
  // - 1 = potentially coalescing in a time intervall or
  // - 2 = potentially recombining in a time intervall
  size_t state_1; // State of active_node_1
  size_t state_2; // -"-      active_node_2

  // Placeholders for the rates at which things happen for the active nodes
  double rate_1;
  double rate_2;
  
  double current_time;
  size_t event_nr;

  // If the start_node is above the local tree, then we start with coalescence
  // of the local root
  if ( start_node->height() > active_node_2->height() ) start_node = active_node_2;

  for (EventIterator event = EventIterator(this, start_node); event.good(); ++event) {
    dout << "* * Time interval: " << (*event).start_height() << " - "
         << (*event).end_height() << std::endl;

    // What is currently happening?
    state_1 = getNodeState(active_node_1, (*event).start_height());
    state_2 = getNodeState(active_node_2, (*event).start_height());
    assert( state_1 != 0 || state_2 != 0 );                   

    // Calc total rate of anything happening in this time interval
    rate_1 = calcRate(active_node_1, state_1, state_2, *event);
    rate_2 = calcRate(active_node_2, state_2, state_1, *event);

    dout << "* * Active Node 1: " << active_node_1 << " State: " << state_1
         << " Rate: " << rate_1 << std::endl;
    dout << "* * Active Node 2: " << active_node_1 << " State: " << state_2 
         << " Rate: " << rate_2 << std::endl;
    assert( rate_1 + rate_2 > 0 );

    // Sample the time at which the next thing happens
    current_time = sampleExpTime(rate_1 + rate_2, (*event).length());
    dout << "* * Next event at time " << current_time
         << " (-1 == no event)" << std::endl;

    // Go on if nothing happens in this time interval
    if ( current_time == -1 ) {
      if (state_1 == 2) active_node_1 = possiblyMoveUpwards(active_node_1, *event);
      if (state_2 == 2) active_node_2 = possiblyMoveUpwards(active_node_2, *event);
      continue;
    }
    
    // Determine for which active node the event occurred
    event_nr = sampleWhichRateRang(rate_1, rate_2);
    dout << "* * Event for active node " << event_nr << std::endl;

    // Pw. coalescence
    if (state_1 == 1 && state_2 == 1) {
      if ( (*event).contemporaries().size() == 0 ||
           random_generator()->sample() * (1 + 2 * (*event).contemporaries().size() ) <= 1 ) {
        
        implementPwCoalescence(active_node_1, active_node_2, current_time);
      }
    }
  }  
}

// Function to determine the state of a (branch above an) local node
// - 0 = off
// - 1 = potentially coalescing
// - 2 = potentially recombining
size_t Forest::getNodeState(Node const *node, const double &current_time) const {
  if (node->height() > current_time) return(0);
  if (node->is_root()) return(1);
  if (!node->local()) return(2);
  return(0);
}

void Forest::implementPwCoalescence(Node* root_1, Node* root_2, const double &time) {
  dout << "* * Both nodes coalesced together" << std::endl;
  dout << "* * Implementing...";
  Node* new_root = NULL;

  // both nodes may or may not mark the end of a single branch at the top of their tree,
  // which we don't need anymore.
  if (root_1->numberOfChildren() == 1) {
    if (root_2->numberOfChildren() == 1) {
      dout << " root_2 has" << std::endl;
      // both trees have a single branch => delete one
      root_2 = root_2->lower_child();
      this->nodes()->remove(root_2->parent());
      assert( root_2 != NULL );
    }
    // (now) only root_1 has a single branch => use as new root
      dout << " root_1 has" << std::endl;
    this->nodes()->move(root_1, time);
    new_root = root_1;
    root_1 = root_1->lower_child();
    assert( root_1 != NULL );
  } 
  else if (root_2->numberOfChildren() == 1) {
      dout << " root_2 has" << std::endl;
    // only root_2 has a single branch => use as new root
    new_root = root_2;
    root_2 = root_2->lower_child();
  }
  else {
      dout << " root_2 has" << std::endl;
    // No tree a has single branch on top => create a new root
    new_root = new Node(time);
    this->nodes()->add(new_root);
  }
  dout << " done" << std::endl;

  root_1->set_parent(new_root);
  root_2->set_parent(new_root);
  new_root->set_higher_child(root_1);
  new_root->set_lower_child(root_2);
  dout << " done" << std::endl;
  new_root->sort_children();
  dout << " done" << std::endl;
}

void Forest::sampleCoalescences2(Node *start_node, const bool &for_initial_tree) {
  assert(start_node->is_root());

  Node* coalescence_root_1 = start_node;
  Node* coalescence_root_2 = this->local_root();
  bool  coalescence_1_local = false;
  bool  coalescence_2_local = false;

  double current_time = std::min(start_node->height(), coalescence_root_2->height());
  double std_expo_sample = -1;
  double expo_sample = -1;

  while (true) {
    EventIterator events = EventIterator(this, coalescence_root_1);
    Event current_event = ++events;
    std_expo_sample = this->random_generator()->sampleExpo(1);

    while (true) {
      assert(this->checkTree());
      dout << "* * new event (time " << current_time << ")" << std::endl;

      coalescence_1_local = current_time >= coalescence_root_1->height();
      coalescence_2_local = current_time >= coalescence_root_2->height();
      assert(coalescence_1_local || coalescence_2_local);

      // If root_1 reaches the height of the local_root, we need to move root_2
      // upwards though the (possibly existing) primary tree,  
      if (coalescence_2_local && coalescence_root_2 == this->local_root()) {
        dout << "* * Passing local_root. Moving Root 2 upwards." << std::endl;
        coalescence_root_2 = moveUpwardsInTree(coalescence_root_2);
        dout << "* * * New Root 2: " << coalescence_root_2 << std::endl;
        coalescence_2_local = current_time >= coalescence_root_2->height();
      }

      dout << "* * Active: ";
      if (coalescence_1_local) dout << coalescence_root_1 << " ";
      if (coalescence_2_local) dout << coalescence_root_2 << " ";
      dout << std::endl;
      

      // If SMC or initial tree:
      // Remove branch above recombination point to avoid back coalescence.
      if ( for_initial_tree  || this->model().is_smc_model() ) {
        dout << "* * removing branch above coalescing node" << std::endl;
        current_event.removeFromContemporaries(start_node);
      }

      // Calculate rate of any coalescence occurring
      double rate = 0.05;
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
      current_event = ++events;
    }
        
    dout << "* * coalescence at time " << current_time << std::endl;
    dout << "* * #contemporaries: " << current_event.contemporaries().size() << std::endl;
 
    Node *coalescing_root = NULL, *coalescence_target = NULL;
    
    if (coalescence_1_local && coalescence_2_local) {
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
      // One local coalescence
      assert(coalescence_1_local || coalescence_2_local);
      if (coalescence_1_local) coalescing_root = coalescence_root_1;
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

    // Root 1 and 2 coalesced together? => we are done
    if (coalescing_root == coalescence_root_1 && coalescence_target == coalescence_root_2) {
      dout << "* * Both roots coalesced together" << std::endl;
      updateAbove(local_root(), false, false);    // To local the old local root
      updateAbove(coalescing_root, false, false); // Update the new node
      break;
    }

    // Root 2 coalesced? Move upwards until we hid another root
    if (coalescing_root == coalescence_root_2) {
      dout << "* * Root 2 coalesced" << std::endl;
      coalescence_root_2 = moveUpwardsInTree(coalescence_root_2);
      coalescence_root_2->make_local();
      dout << "* * New Root 2: " << coalescence_root_2 << std::endl;
    }

    // Root 1 coalesced? Follow branch up until
    // - we hit a local branch or the local root => Done  
    //   (happens iff we originally hit the primary tree below the local root
    // - we hit the primary root  => Done 
    //   (happens iff we hit the primary tree elsewhere
    // - we hit another root      => continue coalescence
    //   (happens iff we hit another tree
    if (coalescing_root == coalescence_root_1) {
      dout << "* * Root 1 coalesced" << std::endl;
      coalescence_root_1 = moveUpwardsInTree(coalescence_root_1);

      // If it is not root, then we hit a local node
      if ( (!coalescence_root_1->is_root()) || coalescence_root_1 == local_root() ) {
        dout << "* * * is local (or local root)" << std::endl;
        break;
      }

      // Else also local the current node and continue coalescing
      dout << "* * * is another root" << std::endl;
      coalescence_root_1->make_local();
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

Node* Forest::possiblyMoveUpwards(Node* node, const Event &event) const {
  if ( node->parent_height() == event.end_height() ) {
    node->set_last_update(this->current_base());
    return node->parent();
  }
  return node;
}

// Starting from "node", this moves upwards in the tree until it hits a
// root or an local node. Updates the nodes on its way up.
// If it encounters a non-local branch, it also checks for possible postponed
// recombination events. If there is one, it cuts the subtree out and returns
// the 
// The final root or local node is returned.
Node* Forest::moveUpwardsInTree(Node* node) {
  while ( !node->is_root() ) {
    updateAbove(node, false, false); // Updates only "node"
    node = node->parent();
    dout << "* * * Moving upwards: " << node << std::endl;        
    if (node->local()) {
      updateAbove(node, false); // Updates the rest of the tree
      return node;
    }
  }

  updateAbove(node, false, false); // Updates the found root (only)
  return node;
}


// Calculates the rate of any coalescence occuring in an intervall with a total of 
// lines_number lines out of which coal_lines_number are not coalescenced.
double Forest::calcRate(Node* node, const int &state, const int &other_state, const Event &event) const {
  // Node is off
  if (state == 0) return 0;

  // Coalescence
  if (state == 1 && other_state != 1)
    return ( event.contemporaries().size() / ( 2.0 * this->model().population_size() ) );
  if (state == 1 && other_state == 1)
    return ( 2 * event.contemporaries().size() + 0.5 ) / ( 2.0 * this->model().population_size() );

  // Recombination
  if (state == 2)
    return ( model().recombination_rate() * (this->current_base() - node->last_update()) );

  throw std::logic_error("Error calculating rates");
}

// Looks whether an exp(rate) distributed waiting time runs off in a time
// interval of a given length. Return the time if so, and -1 otherwise.
double Forest::sampleExpTime(double rate, double intervall_length) {
  if (this->expo_sample_ == -1) expo_sample_ = this->random_generator()->sampleExpo(1);

  // Scale interval with rate to bring it on exp(1) timescale
  intervall_length *= rate;
  if (intervall_length < expo_sample_) {
    expo_sample_ -= intervall_length; // Still exp(1) distributed
    return -1;
  } else {
    double time = expo_sample_;
    expo_sample_ = -1;
    return (time/rate);
  }
}
  
size_t Forest::sampleWhichRateRang(const double &rate_1, const double &rate_2) const {
  if (rate_1 == 0) return 2;
  if (rate_2 == 0) return 1;

  if (this->random_generator()->sample() * (rate_1 + rate_2) > rate_1) return 2;
  return 1;
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
  ultimate_root->make_nonlocal(current_base());
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

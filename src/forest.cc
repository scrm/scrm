#include "forest.h"

/******************************************************************
 * Constructors & Initialization
 *****************************************************************/

Forest::Forest() {
  this->initialize();
};

Forest::Forest(Model* model, RandomGenerator* random_generator) {
  this->initialize(model, random_generator);
}

// Sets member variable to default values
void Forest::initialize(Model* model, 
                        RandomGenerator* rg) {

  this->nodes_ = new NodeContainer();
  this->set_model(model);
  this->set_random_generator(rg);
  this->set_current_base(1);
  this->expo_sample_ = -1;
  this->prune_countdown_ = model->prune_interval();
}



/******************************************************************
 * Destructor
 *****************************************************************/

Forest::~Forest() { 
  dout<<"Forest destructor is called"<<endl;
  //delete[] nodes_;
  nodes()->clear();
  dout<<"Forest is deleted"<<endl;
}


/** 
 * function that cuts a subtree out of a tree of the forest and reinserts it as
 * a separate tree.
 *
 * This is primarily used to cut the subtree below an recombination
 * away. 
 *
 * \param cut_point   A TreePoint marking the top of the subtree to cut.
 * \return            The root of the now separated subtree.
 */
Node* Forest::cut(const TreePoint &cut_point) {
  //The node above the cut_point in the old tree
  Node* parent = cut_point.base_node()->parent();
  assert( parent != NULL );

  //The new end of the old branch after the cut
  Node* new_leaf = new Node(cut_point.height(), false, current_base(), 0, 0);
  if ( !cut_point.base_node()->local() )
    new_leaf->set_last_update(cut_point.base_node()->last_update() );
  new_leaf->set_parent(parent);
  parent->change_child(cut_point.base_node(), new_leaf);
  nodes()->add(new_leaf, cut_point.base_node());
  // Update all local nodes above.
  updateAbove(parent, false, true, true);
  dout << "* * New leaf of local tree: " << new_leaf << std::endl;

  //The new "root" of the newly formed tree
  Node* new_root = new Node(cut_point.height(), 
                            cut_point.base_node()->local());
  cut_point.base_node()->set_parent(new_root);
  new_root->set_first_child(cut_point.base_node());
  nodes()->add(new_root, new_leaf);
  dout << "* * New root of subtree: " << new_root << std::endl;

  dout << "* * Done" << std::endl;
  return(new_root);
}


/**
 * Function to update the invariants (local, samples_below, length_below) 
 * of a 'node' and all of its (grand-)parents. Also sets local_root_ if it
 * encounters it. 
 *
 *  \param node       The node at which the functions starts updating the
 *                    invariants. Then updates it's parent and the parents
 *                    parent.
 *  \param above_local_root If true, it uses a faster algorithm that is only correct
 *                    for nodes above the local root. Default false. Best don't touch
 *                    this.
 *  \param recursive  If false, only the given node is updated, but not its parent.
 *                    Default true.
 *  \param dont_localize If true, the function does not set non-local nodes local.
 *                    Needed for updating the tree above the 'left_over' above
 *                    a recombination.
 */
// Think this function should not make any non-local local, because we may
// ignore postponed recombinations this way. Will investigate this further. 
// -Paul
void Forest::updateAbove(Node* node, bool above_local_root, bool recursive, bool dont_localize) {
  //dout << "* * * Updating: " << node << " local: " << node->local()
  //     << " fastforward: " << above_local_root << std::endl;

  // Fast forward above local root because this part is non-local
  if (above_local_root) {
    if (node->local()) node->make_nonlocal(current_base());
    node->set_samples_below(this->sample_size());
    node->set_length_below(this->local_root()->length_below());
    if ( node->is_root() ) {
      set_primary_root(node);
      return;
    }
    if ( recursive ) updateAbove(node->parent(), true);
    return;
  }

  // Calculate new values for samples_below and length_below for the current
  // node
  Node *l_child = node->first_child();
  Node *h_child = node->second_child();

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

  // If nothing changed, we also don't need to update the tree further above...
  // (Is that really the case?)
  if (recursive &&
      samples_below == node->samples_below() && 
      areSame(length_below, node->length_below()) ) {
    //dout << "FF STOPED:" << node << std::endl; 
    return;
  }

  // Update the node
  if ( samples_below == 0 ) node->make_nonlocal(current_base());
  else if ( samples_below < this->sample_size() && !dont_localize ) node->make_local();
  node->set_samples_below(samples_below);
  node->set_length_below(length_below);

  // Check if we are above or equal the local root
  if ( samples_below == this->sample_size() ) {
    node->make_nonlocal(current_base());
    // Are we the local root?
    if (l_child->samples_below() > 0 && h_child->samples_below() > 0) {
      //dout << "* * * is local_root" << std::endl;
      set_local_root(node);
    }
    if ( node->is_root() ) set_primary_root(node);
    above_local_root = true;
  }

  // Go further up if possible
  if ( recursive && !node->is_root() ) {
    updateAbove(node->parent(), above_local_root, recursive, dont_localize);
  }
}


/**
 * Function that build the initial tree at the very left end of the sequence.
 *
 * Also creates the sample nodes.
 */
void Forest::buildInitialTree() {
  dout << "===== BUILDING INITIAL TREE =====" << std::endl;
  dout << "* Adding first node... ";
  Node* first_node = new Node(0, true, 0, 1);
  this->nodes()->add(first_node);
  this->set_local_root(first_node);
  this->set_primary_root(first_node);
  first_node->set_label(1);
  dout << "done." << std::endl;

  for (size_t i=1; i < this->model().sample_size(); i++) {
    this->set_sample_size(i+1);
    dout << "* adding node ";
    //Create a new separate little tree of and at height zero
    Node* new_leaf = new Node(0, true, 0, 1);
    new_leaf->set_label(i+1);
    dout << "(" << new_leaf << ")" << std::endl;
    nodes()->add(new_leaf);
    dout << "* staring coalesces" << std::endl;

    //Coalesces the separate tree into the main tree
    this->sampleCoalescences(new_leaf, false);
    dout << "* * Tree:" << std::endl;
    assert(this->printTree());
    assert(this->printNodes());
    assert(this->checkTree());
  }
  this->set_next_base();
  writeTree(this->local_root(),this->model_->population_size(), 0.0);
  //cout<<this->local_root()->tree_topo_bl<<endl;
  
    //set the index for all forest nodes....
  for(size_t i = 0; i < this->getNodes()->size(); ++i) {
	this->nodes()->get_copy(i)->index = i;  
  }
  //size_t i = 0;
  //for (ConstNodeIterator it = this->nodes()->iterator(); it.good(); ++it) {this->nodes()->get_copy(it)->index = i; i++;};
  //Node * current_node=this->getNodes()->first();
  //size_t i = 0;
  //while (i!= this->getNodes()->size()-1){
	   //current_node->index=i;
	  //current_node = current_node->next();
	//i++;
  //}
  	//current_node->index=i;
	
  
}


/**
 * Uniformly samples a TreePoint on the local tree.
 *
 * Its arguments are meant to be used only when the function iteratively calls
 * itself. Just call it without any arguments if you want to sample a TreePoint.
 *
 * The function first samples a part of the total height of the tree and then
 * goes down from the root, deciding at each node if that point is to the left
 * or right, which should give us an O(log(#nodes)) algorithm.
 *
 * I checked the distribution of this function in multiple cases. -Paul
 *
 * \param node The current position in the tree when the functions goes down
 *             iteratively.
 *
 * \param length_left The length that is left until we encounter the sampled
 *              length.
 *
 * \return The sampled point on the tree.
 */
TreePoint Forest::samplePoint(Node* node, double length_left) {
  if (node == NULL) {
    // Called without arguments => initialization
    node = this->local_root();
    length_left = this->random_generator()->sample() * node->length_below();
  }

  assert( node->local() || node == this->local_root() );
  assert( length_left >= 0 );
  assert( length_left < (node->length_below() + node->height_above()) );

  if ( node != this->local_root() ) {
    if ( length_left < node->height_above() ) {
      return TreePoint(node, length_left, true);
    }

    length_left = length_left - node->height_above();
    assert( length_left >= 0 );
  }

  // At this point, we should have at least one local child
  assert( node->first_child() != NULL );
  assert( node->first_child()->local() || node->second_child()->local() );

  // If we have only one local child, then give it the full length we have left.
  if ( !node->first_child()->local() ) {
    return samplePoint(node->second_child(), length_left);
  }
  if ( node->second_child() == NULL || !node->second_child()->local() ) {
    return samplePoint(node->first_child(), length_left);
  }

  // If we have two local children, the look if we should go down left or right.
  double tmp = node->first_child()->height_above() + node->first_child()->length_below();
  if ( length_left <= tmp )
    return samplePoint(node->first_child(), length_left);
  else 
    return samplePoint(node->second_child(), length_left - tmp);
}


/** 
 * Function to modify the tree after we encountered a recombination on the
 * sequence. Also samples a place for this recombination on the tree, marks the
 * branch above as non-local (and updates invariants) if needed, cuts the
 * subtree below away and starts a coalescence from it's root.
 */
void Forest::sampleNextGenealogy() {
  this->set_current_base(next_base_);
  dout << std::endl << "===== BUILDING NEXT GENEALOGY =====" << std::endl;
  dout << "Sequence position: " << this->current_base() << std::endl;

  // Check if we prune while building the genealogy
  pruning_ = false;
  if ( model().exact_window_length() != 0 && current_base() >= model().exact_window_length() ) {
    if (prune_countdown_ == 0) {
      prune_countdown_ = model().prune_interval();
      pruning_ = true;
    }
    --prune_countdown_;
  }
  if (pruning_) dout << "We will prune the tree this time" << std::endl;

  // Sample the recombination point
  TreePoint rec_point = this->samplePoint();

  dout << "* Recombination at height " << rec_point.height() << " ";
  dout << "(above " << rec_point.base_node() << ")"<< std::endl;

  dout << "* Cutting subtree below recombination " << std::endl;
  this->cut(rec_point);

  assert( rec_point.base_node()->local() );
  assert(this->printTree());

  dout << "* Starting coalescence" << std::endl;
  this->sampleCoalescences(rec_point.base_node()->parent(), pruning_);



  ////set the index for all forest nodes....
  for(size_t i = 0; i < this->getNodes()->size(); ++i) {
	this->nodes()->get_copy(i)->index = i;  
  }
//Node * current_node=this->getNodes()->first();
  //size_t i = 0;
  //while (i!= this->getNodes()->size()-1){
	   //current_node->index=i;
	  //current_node = current_node->next();
	//i++;
  //}
  	//current_node->index=i;
  	

  assert(this->printTree());
  assert(this->printNodes());
  assert(this->checkTree());
  this->set_next_base();
  writeTree(this->local_root(),this->model_->population_size(),0);
}


/** 
 * Function for doing a coalescence.
 *
 * \param start_node The node at which the coalescence starts. Must be the root
 *                   of a tree.
 */
void Forest::sampleCoalescences(Node *start_node, bool pruning) {
  assert( start_node->is_root() );
  // We can have one or active local nodes: If the coalescing node passes the
  // local root, it also starts a coalescence.
  Node* active_node_1 = start_node;
  Node* active_node_2 = this->local_root();

  // States: Each (branch above a) node can either be in state
  // - 0 = off (the other coalescence has not reached it yet) or
  // - 1 = potentially coalescing in a time interval or
  // - 2 = potentially recombining in a time interval
  size_t state_1; // State of active_node_1
  size_t state_2; // -"-      active_node_2

  // Placeholders for the rates at which things happen for the active nodes
  double rate_1;
  double rate_2;

  double current_time = -1;
  size_t event_nr;
  Node **event_node, **other_node;
  size_t event_state, other_state;
  TreePoint event_point;

  // If the start_node is above the local tree, then we start with coalescence
  // of the local root
  if ( start_node->height() > active_node_2->height() ) start_node = active_node_2;

  for (TimeIntervalIterator event = TimeIntervalIterator(this, start_node, pruning); 
       event.good(); ++event) {
    
    dout << "* * Time interval: " << (*event).start_height() << " - "
        << (*event).end_height() << std::endl;

    assert( current_time < 0 || current_time == (*event).start_height() );

    // What is currently happening?
    state_1 = getNodeState(active_node_1, (*event).start_height());
    state_2 = getNodeState(active_node_2, (*event).start_height());

    // Calc total rate of anything happening in this time interval
    rate_1 = calcRate(active_node_1, state_1, state_2, *event);
    rate_2 = calcRate(active_node_2, state_2, state_1, *event);

    dout << "* * * Active Node 1: " << active_node_1 << " State: " << state_1
        << " Rate: " << rate_1 << std::endl;
    dout << "* * * Active Node 2: " << active_node_2 << " State: " << state_2 
        << " Rate: " << rate_2 << std::endl;

    assert( active_node_1 != active_node_2 );
    assert( state_1 != 0 || state_2 != 0 );
    assert( state_1 == 1 || active_node_1->parent_height() >= current_time );
    assert( state_2 == 1 || active_node_2->parent_height() >= current_time );

    // Sample the time at which the next thing happens
    current_time = sampleExpTime(rate_1 + rate_2, (*event).length());
    if (current_time >= 0) {
      current_time += (*event).start_height();
      dout << "* * Event at time " << current_time << std::endl;
    }

    // Go on if nothing happens in this time interval
    if ( current_time == -1 ) {
      if (state_1 == 2) {
        active_node_1 = possiblyMoveUpwards(active_node_1, *event);
        if (active_node_1->local()) {
          dout << "* * * Active Node 1 hit a local node. Done" << std::endl;
          updateAbove(active_node_1);
          return;
        }
      }
      if (state_2 == 2) active_node_2 = possiblyMoveUpwards(active_node_2, *event);

      if (active_node_1 == active_node_2) {
        dout << "* * * Active Nodes hit each other in " << active_node_1 
            << ". Done." << std::endl;
        updateAbove(active_node_1);
        return;
      }
      continue;
    }

    // First take care of pairwise coalescence
    if (state_1 == 1 && state_2 == 1) {
      if ( (*event).numberOfContemporaries() == 0 ||
          random_generator()->sample() * (1 + 2 * (*event).numberOfContemporaries() ) <= 1 ) {

        implementPwCoalescence(active_node_1, active_node_2, current_time);
        return;
      }
    }

    // Now look for coalescence of only one line and recombinations
    //
    // First determine for which active node the event occurred, so we
    // don't have to duplicate the code.
    event_nr = sampleWhichRateRang(rate_1, rate_2);
    if (event_nr == 1) {
      event_node = &active_node_1;
      event_state = state_1;
      other_node = &active_node_2;
      other_state = state_2; 
    } else {
      event_node = &active_node_2;
      event_state = state_2;
      other_node = &active_node_1;
      other_state = state_1; 
    }

    // Now do the real work
    assert( event_state != 0 );
    if ( event_state == 1 ) {
      // Coalescence: sample target point and implement the coalescence
      dout << "* * * Active Node " << event_nr  << ": Coalescence" << std::endl;
      dout << "* * * #Contemporaries: " << (*event).numberOfContemporaries() << std::endl;
      event_point = TreePoint((*event).getRandomContemporary(), current_time, false);
      dout << "* * * Above node " << event_point.base_node() << std::endl;
      *event_node = implementCoalescence(*event_node, event_point);

      // If the other node was looking for a recombination, we must ensure that
      // the branch below the event get marked as updated.
      if ( other_state == 2 ) {
        // If the coalescing node coalesced into the branch directly above 
        // the recombining node, then we are done.
        if ( (*other_node)->parent() == *event_node ) {
          (*other_node)->set_last_update(this->current_base());
          dout << "* * * Recombining Node moved into coalesced node. Done." << std::endl;
          updateAbove(*other_node, false, false);
          updateAbove((*other_node)->parent());
          return;
        }
      }

      // Check if are can stop.
      if ( (*event_node)->local() ) {
        // Only one node should be active if there are still local nodes at the
        // current height. Hence we are done if so.
        assert( other_state == 0 );

        dout << "* * * We hit the local tree. Done." << std::endl;
        updateAbove(*event_node); 
        return;
      }

      if (active_node_1 == active_node_2) {
        dout << "* * * Coalescend into other Active Node. Done." << std::endl;
        updateAbove(*event_node); 
        return;
      }

      // If we hit an non-local branch:
      // Begin next interval at the coalescence height and remove the branch
      // below from contemporaries.
      event.splitCurrentInterval(*event_node, event_point.base_node());
      assert(this->printTree());
    }
    else if (event_state == 2) {
      // Recombination: sample point of recombination and implement
      dout << "* * * Active Node " << event_nr<< ": Recombination" << std::endl;
      updateAbove(*event_node, false, false);
      event_point = TreePoint(*event_node, current_time, false);
      *event_node = cut(event_point);
      updateAbove(*event_node, false, false);
      event.recalculateInterval();

      assert(this->printTree());
      continue;
    }
  }  
}


/**
 *  Function to determine the state of (the branch above) a node for an ongoing
 *  coalescence.
 *
 *  The States are: 0 = off, 1 = potentially coalescing, 2 = potentially recombining
 *
 *  \param node the node for which to tell the start
 *  \param current_time the time at which the coalescence is
 *  \return The state of the node
 */
size_t Forest::getNodeState(Node const *node, const double &current_time) const {
  if (node->height() > current_time) return(0);
  if (node->is_root()) return(1);
  if (!node->local()) return(2);
  dout << "Error getting node state." << std::endl;
  dout << "Height: " << node->height() 
       << " current time: " << current_time 
       << " diff: " << node->height() - current_time << std::endl;
  dout << "Node local: " << node->local() << std::endl;
  dout << "Node root: " << node->is_root() << std::endl;
  assert( false );
  return(0);
}


/**
 * Modifies the forest to reflect the coalescences of a coalescing node into
 * a point on a existing tree.
 *
 * Attention: The returned node does still require an update!
 *
 * \param coal_node The coalescing node
 * \param coal_point The point on the tree into which coal_node coalesces.
 *
 * \return The new node at the point of coalescence.
 */
Node* Forest::implementCoalescence(Node* coal_node, const TreePoint &coal_point) {
  Node* new_node;

  // Look if we can reuse the root that coalesced as new node
  if ( coal_node->numberOfChildren() == 1 ){
    new_node = coal_node;
    coal_node = coal_node->first_child();
    nodes()->move(new_node, coal_point.height());
  } else {
    new_node = new Node(coal_point.height());
    new_node->change_child(NULL, coal_node);
    coal_node->set_parent(new_node);
    nodes()->add(new_node);
  }

  // Now we have:
  // coal_point.base_node() = Node in the target tree under the coalescence
  // coal_node =  Root of the coalescing tree 
  // new_node =   New parent of both other nodes

  // Update new_node
  new_node->change_child(NULL, coal_point.base_node());
  new_node->set_parent(coal_point.base_node()->parent());
  if (!coal_point.base_node()->local()) {
    new_node->set_local(false);
    new_node->set_last_update(coal_point.base_node()->last_update());
  } else {
    new_node->set_local(true);
  }

  // Integrate it into the tree
  coal_point.base_node()->set_parent(new_node);
  new_node->parent()->change_child(coal_point.base_node(), new_node);

  // And update
  updateAbove(coal_node, false, false); // Update coal_node ONLY
  // updateAbove(new_node, false, false); 
  // Updating the new_node will always activate it, but there may still be
  // unimplemented recombinations above

  return(new_node);
}


/** 
 * Modifies the forest to reflect that two coalescing nodes coalesced together.
 * 
 * \param root_1 The first coalescing node
 * \param root_2 The second coalescing node
 * \param time   The time at which the coalescence happens
 *
 */
void Forest::implementPwCoalescence(Node* root_1, Node* root_2, const double &time) {
  dout << "* * Both nodes coalesced together" << std::endl;
  dout << "* * Implementing..." << std::flush;
  Node* new_root = NULL;

  // both nodes may or may not mark the end of a single branch at the top of their tree,
  // which we don't need anymore.
  if (root_1->numberOfChildren() == 1) {
    if (root_2->numberOfChildren() == 1) {
      // both trees have a single branch => delete one
      root_2 = root_2->first_child();
      this->nodes()->remove(root_2->parent());
      assert( root_2 != NULL );
    }
    // (now) only root_1 has a single branch => use as new root
    this->nodes()->move(root_1, time);
    new_root = root_1;
    root_1 = root_1->first_child();
    assert( root_1 != NULL );
  } 
  else if (root_2->numberOfChildren() == 1) {
    // only root_2 has a single branch => use as new root
    this->nodes()->move(root_2, time);
    new_root = root_2;
    root_2 = root_2->first_child();
  }
  else {
    // No tree a has single branch on top => create a new root
    new_root = new Node(time);
    this->nodes()->add(new_root);
  }

  root_1->set_parent(new_root);
  root_2->set_parent(new_root);
  new_root->set_second_child(root_1);
  new_root->set_first_child(root_2);

  updateAbove(root_1, false, false);
  updateAbove(root_2, false, false);
  updateAbove(new_root, false, false);
  dout << " done" << std::endl;
}


/** 
 * Helper function for doing a coalescence.
 * Moves the 'active' flag (i.e. the node stored in root_1 or root_2 in sampleCoalescence)
 * from a node to it's parent if the branch above the node
 * ends this the current time interval.
 *
 * This function is used the pass the active flag upwards in the tree if the
 * node is active, but neither coalescing nor a recombination happens on the
 * branch above, e.g. after a local-branch became active because it was hit by a
 * coalescence or a non-local branch was active and no recombination occurred.
 *
 * Also updates the active node if it moves up.
 *
 * \param node An active node
 * \param time_interval The time interval the coalescence is currently in.
 * 
 * \return  Either the parent of 'node' if we need to move upwards or 'node'
 *          itself
 */
Node* Forest::possiblyMoveUpwards(Node* node, const TimeInterval &time_interval) {
  if ( node->parent_height() == time_interval.end_height() ) {
    node->set_last_update(this->current_base());
    updateAbove(node, false, false);
    return node->parent();
  }
  return node;
}


/** 
 * Calculates the rate of an event occurring above a node in an time interval.
 *
 * In case both nodes are coalescing, it splits the rate of pair wise
 * coalescence equally on both nodes.
 *
 * \param node The node for which we calculation the rate.
 * \param state The state of 'node'
 * \param other_state The state of the other active node
 * \param event The current time interval
 *
 * \returns The rate.
 */
double Forest::calcRate(Node* node, 
                        const int &state, 
                        const int &other_state, 
                        const TimeInterval &event) const {
  // Node is off
  if (state == 0) return 0;

  // Coalescence
  if (state == 1 && other_state != 1)
    return ( event.numberOfContemporaries() / ( 2.0 * this->model().population_size() ) );
  if (state == 1 && other_state == 1)
    return ( 2 * event.numberOfContemporaries() + 0.5 ) / ( 2.0 * this->model().population_size() );

  // Recombination
  if (state == 2)
    return ( model().recombination_rate() * (this->current_base() - node->last_update()) );

  throw std::logic_error("Error calculating rates");
}


/**
 * Looks whether an exponentially distributed waiting time runs off in a time
 * interval. 
 *
 * This function reuses the same sample for multiple interval as long as the
 * waiting time does not run off, which greatly improves performance and is
 * correct because of the Markov property of waiting times.
 *
 * \param rate The rate of the waiting time
 * \interval_length The length of the current time interval.
 *
 * \return The time at which the rate ran of if it did, or -1 otherwise.
 */
double Forest::sampleExpTime(double rate, double interval_length) {
  if (rate == 0) return -1;
  if (this->expo_sample_ == -1) expo_sample_ = this->random_generator()->sampleExpo(1);

  // Scale interval with rate to bring it on exp(1) timescale
  interval_length *= rate;
  if (interval_length < expo_sample_) {
    expo_sample_ -= interval_length; // Still exp(1) distributed
    return -1;
  } else {
    double time = expo_sample_;
    expo_sample_ = -1;
    return (time/rate);
  }
}


/**
 * Sample which of two waiting time ran off, conditioned on that one of them
 * did.
 *
 * \param rate_1 The rate of the first waiting time
 * \param rate_2 The rate of the second waiting time
 *
 * \return 1 if the first one ran off, and 2 if the second one did.
 */
size_t Forest::sampleWhichRateRang(const double &rate_1, const double &rate_2) const {
  if (rate_1 == 0) return 2;
  if (rate_2 == 0) return 1;

  if (this->random_generator()->sample() * (rate_1 + rate_2) > rate_1) return 2;
  return 1;
}


/**
 * Test if a branch above a given node can be pruned, i.e. removed from the
 * forest.
 *
 * Currently, there are three types of nodes that will be pruned:
 * 1) orphaned roots: roots with no children, i.e. a tree consisting just of one
 * node. This single node trees are created when all other nodes on the tree are
 * pruned.
 * 2) in-between nodes: nodes that artificially split a single branch into two
 * parts. Are created if the former other child of the node is pruned.
 * 3) old nodes: If specified, non-local nodes that have not been updated for a
 * while can be pruned to keep the tree from exploding in size.
 * 
 * Gets called from the TimeIntervalIterator when creating time intervals.
 *
 * \param node The node to be examined for pruning
 * 
 * \param return TRUE iff the node can be pruned
 */
bool Forest::isPrunable(Node const* node) const {
  if ( model().exact_window_length() == 0 ) return false;
  if ( node == local_root() ) return false;

  switch ( node->numberOfChildren()) {
    case 2: return false;

    case 1:
      // in-between branches (= 1 child, 1 parent (=not root), updated at same time)
      if ( (!node->is_root()) &&
            node->last_update() == node->first_child()->last_update() ) return true;
      return false;

    case 0:
      // remove orphaned roots
      if ( node->is_root() ) return true;  

      // remove old external nodes
      if ( (!node->local()) && 
          this->current_base() - node->last_update() > model().exact_window_length() ) { 
        return true;
      }
      return false;
  }
  assert(0);
  return false; 
}


/**
 * Remove a prunable node from the tree.
 *
 * This function perform the actual pruning. See isPunable for details about
 * pruning.
 *
 * Gets called from the TimeIntervalIterator when creating time intervals.
 *
 * \param node The node to be pruned.
 */
void Forest::prune(Node *node) {
  assert( isPrunable(node) );
  assert( ! (node->is_root() && node->local()) );
  assert( !node->in_sample() );

  dout << "* * * PRUNING: Removing node " << node << " from tree" << std::endl;
  dout << "* * * PRUNING: Parent: " << node->parent() << "; " 
       << "l_child: " << node->first_child() << "; "
       << "h_child: " << node->second_child() << "; "
       << "local: " << node->local() << std::endl;

  /* Possible Cases:
   * + Orphaned root => just delete 
   * + In-between node => relink and delete
   * + Old branch => unset in parent and delete 
   * */

  // In-between node
  if ( node->numberOfChildren() == 1 ) {
    dout << "* * * PRUNING: Reason: Unneeded Node" << std::endl;
    Node* child = node->first_child();
    child->set_parent( node->parent() );
    node->parent()->change_child(node, child); 
  } 

  // External branch 
  else if ( !node->is_root() ) {
    dout << "* * * PRUNING: Reason: Old external branch" << std::endl;
    node->parent()->remove_child(node);
  }

  // And delete
  nodes()->remove(node);
}

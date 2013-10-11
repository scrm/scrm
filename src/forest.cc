/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013 Paul R. Staab, Sha (Joe) Zhu and Gerton Lunter
 * 
 * This file is part of scrm.
 * 
 * scrm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
 
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

  this->set_model(model);
  this->set_random_generator(rg);
  this->set_current_base(1);
  this->prune_countdown_ = model->prune_interval();
  this->set_sample_size(0);
}



/******************************************************************
 * Destructor
 *****************************************************************/

Forest::~Forest() { 
  //dout<<"Forest destructor is called"<<endl;
  nodes()->clear();
  //dout<<"Forest is deleted"<<endl;
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
  new_leaf->set_population(cut_point.base_node()->population());
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
  new_root->set_population(cut_point.base_node()->population());
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
 * Function that builds the initial tree at the very left end of the sequence.
 *
 * Also creates the sample nodes.
 */
void Forest::buildInitialTree() {
  dout << "===== BUILDING INITIAL TREE =====" << std::endl;
  dout << "* Adding first node... ";
  Node* first_node = new Node(model().sample_time(1), true, 0, 1, 0, 1);
  first_node->set_population(model().sample_population(1));
  this->nodes()->add(first_node);
  this->set_local_root(first_node);
  this->set_primary_root(first_node);
  dout << "done." << std::endl;

  for (size_t i=1; i < this->model().sample_size(); i++) {
    this->set_sample_size(i+1);
    
    dout << "* adding node ";
    //Create a new separate little tree of and at height zero
    Node* new_leaf = new Node(model().sample_time(i), true, 0, 1, 0, i+1);
    new_leaf->set_population(model().sample_population(i));
    dout << "(" << new_leaf << ")" << std::endl;
    nodes()->add(new_leaf);
    dout << "* starting coalescences" << std::endl;

    //Coalesces the separate tree into the main tree
    this->sampleCoalescences(new_leaf, false);
    //this->clear_initial_coalevent();
    dout << "* * Tree:" << std::endl;

    assert(this->printNodes());
    assert(this->checkTree());
    assert(this->checkLeafsOnLocalTree());
    assert(this->printTree());
  }
  this->set_next_base();
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
      assert( node->local() );
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
  if ( model().exact_window_length() != -1 && current_base() >= model().exact_window_length() ) {
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
  assert( this->printTree() );

  this->initialize_recomb_coalescent(rec_point.height());
  dout << "* Starting coalescence" << std::endl;
  this->sampleCoalescences(rec_point.base_node()->parent(), pruning_);

  assert( this->printTree() );
  assert( this->checkLeafsOnLocalTree() );
  assert( this->printNodes() );
  assert( this->checkTree() );

  this->set_next_base();
  //writeTree(this->local_root(),this->model_->population_size(),0);
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
  set_active_node(0, start_node);
  set_active_node(1, this->local_root());

  // Placeholders for the rates at which things happen for the active nodes

  tmp_event_ = Event(-1);
  coalescence_finished_ = false;
  this->initialize_event(start_node->height());
  // If the start_node is above the local tree, then we start with coalescence
  // of the local root
  if ( start_node->height() > active_node(1)->height() ) start_node = active_node(1);

  for (TimeIntervalIterator ti(this, start_node, pruning); ti.good(); ++ti) {
    dout << "* * Time interval: " << (*ti).start_height() << " - "
        << (*ti).end_height() << std::endl;

    assert( tmp_event_.time() < 0 || tmp_event_.time() == (*ti).start_height() );

    // Update States & Rates (see their declaration for explanation); 
    states_[0] = getNodeState(active_node(0), (*ti).start_height());
    states_[1] = getNodeState(active_node(1), (*ti).start_height());

    // Fixed time events (e.g pop splits/merges & single migration events first
    if (model().hasFixedTimeEvent((*ti).start_height())) implementFixedTimeEvent(ti);

    calcRates(*ti);

    dout << "* * * Active Nodes: a0:" << active_node(0) << ":s" << states_[0]
        << "(p" << active_node(0)->population() << ")" 
        << " a1:" << active_node(1) << ":s" << states_[1] 
        << "(p" << active_node(1)->population() << ")" << std::endl
        << "* * * Total Rates: " << rates_[0] << " " 
        << rates_[1] << " " << rates_[2] << std::endl;

    assert( active_node(0) != active_node(1) );
    assert( states_[0] != 0 || states_[1] != 0 );
    assert( states_[0] == 1 || active_node(0)->parent_height() >= tmp_event_.time() );
    assert( states_[1] == 1 || active_node(1)->parent_height() >= tmp_event_.time() );
    assert( checkContemporaries(*ti) );

    // Sample the time at which the next thing happens
    
    sampleEvent(*ti, tmp_event_time_, tmp_event_line_, tmp_event_);
    dout << "* * * " << tmp_event_ << std::endl;

    // Go on if nothing happens in this time interval
    if ( tmp_event_.isNoEvent() ) {
      if (states_[0] == 2) {
        set_active_node(0, possiblyMoveUpwards(active_node(0), *ti));
        if (active_node(0)->local()) {
          dout << "* * * Active Node 0 hit a local node. Done" << std::endl;
          updateAbove(active_node(0));
          return;
        }
      }
      // There are no local node above the local root, which is the lowest node
      // that active_node(1) can be.
      if (states_[1] == 2) set_active_node(1, possiblyMoveUpwards(active_node(1), *ti));

      if (active_node(0) == active_node(1)) {
        dout << "* * * Active Nodes hit each other in " << active_node(0) 
            << ". Done." << std::endl;
        updateAbove(active_node(0));
        return;
      }
      continue;
    }

    // First take care of pairwise coalescence
    else if ( tmp_event_.isPwCoalescence() ) {
	// Record this  interval (coalescent)
   	   this->record_coalevent((*ti),tmp_event_.time());
	   //this->initialize_event((*ti), tmp_event_.time());    
	   //this->record_coalevent();
        implementPwCoalescence(active_node(0), active_node(1), tmp_event_.time());
      return;
    }

    else if ( tmp_event_.isRecombination() ) {
		// Record this  interval (recombination)
	   this->record_recombevent((*ti),tmp_event_.time());
	   //this->initialize_event((*ti), tmp_event_.time());    
	   //this->record_recombevent();
      this->implementRecombination(tmp_event_, ti);
      continue;
    }

    else if ( tmp_event_.isMigration() ) {
      this->implementMigration(tmp_event_, ti);
      continue;
    }

    else if ( tmp_event_.isCoalescence() ) {
      	// Record this  interval (coalescent)
	   this->record_coalevent((*ti),tmp_event_.time());
	   //this->initialize_event((*ti), tmp_event_.time());    
	   //this->record_coalevent();
      this->implementCoalescence(tmp_event_, ti);
      if (coalescence_finished_) return;

      assert( this->printTree() );
    }
  }
}  


void Forest::calcRates(const TimeInterval &ti) {
  rates_[0] = 0.0;
  rates_[1] = 0.0;
  rates_[2] = 0.0;
  active_nodes_timelines_[0] = 0;
  active_nodes_timelines_[1] = 0;

  // Set rate of first node
  if (states_[0] == 1) {
    rates_[0] += model().total_migration_rate(active_node(0)->population());
    if (model().growth_rate(active_node(0)->population()) == 0.0) 
      rates_[0] += calcCoalescenceRate(active_node(0)->population(), ti); 
    else {
      rates_[1] += calcCoalescenceRate(active_node(0)->population(), ti);
      active_nodes_timelines_[0] = 1;
    }
  }
  else if (states_[0] == 2) {
    rates_[0] += calcRecombinationRate(active_node(0));
  }

  // The second node is a bit more complicated 
  if (states_[1] == 1) {
    rates_[0] += model().total_migration_rate(active_node(1)->population());
    if (model().growth_rate(active_node(1)->population()) == 0.0) {
      // No Growth => Normal time
      rates_[0] += calcCoalescenceRate(active_node(1)->population(), ti);

      if (states_[0] == 1 && active_node(0)->population() == active_node(1)->population()) { 
        // Also add rates for pw coalescence
        rates_[0] += calcPwCoalescenceRate(active_node(1)->population()); 
      }
    }
    else {
      // Growth => we need a exponential time
      if (states_[0] == 1 && active_node(0)->population() == active_node(1)->population()) {
        // We can use the timeline of the first node
        rates_[1] += calcCoalescenceRate(active_node(1)->population(), ti);
        // And must add pw coalescence again
        rates_[1] += calcPwCoalescenceRate(active_node(1)->population()); 
        active_nodes_timelines_[1] = 1;
      } 
      else {
        // We need our own timeline (We are ignoring that both populations could
        // have equal growth rates at the moment...)
        rates_[2] += calcCoalescenceRate(active_node(1)->population(), ti);
        active_nodes_timelines_[1] = 2;
      }
    }
  }
  else if (states_[1] == 2) {
    rates_[0] += calcRecombinationRate(active_node(1));
  }
}  


/**
* Samples if an event (coalescence, recombination or migration of active nodes)
* happens in the current TimeInterval or not. 
*
* In particular requires that the 'temporary' forest members samples_, rates_
* and active_nodes_ are set correctly beforehand.  
*
* \param ti The current time interval
* \returns the event that has happened (can also be a "NoEvent" event)
*/
void Forest::sampleEvent(const TimeInterval &ti, double tmp_event_time, 
                         size_t tmp_event_line, Event &return_event) const {
  tmp_event_time = -1;
  tmp_event_line = -1;

  // Sample on which time and time line the event happens (if any)
  for (size_t i = 0; i < 3; ++i) {
    if (rates_[i] == 0.0) continue;
    selectFirstTime(random_generator()->sampleExpoExpoLimit(rates_[i], getTimeLineGrowth(i), ti.length()), 
                    i, tmp_event_time, tmp_event_line );
  }

  // Correct the time from relative to the time interval to absolute 
  if (tmp_event_time != -1) tmp_event_time += ti.start_height();
  assert( (tmp_event_time == -1) || 
          (ti.start_height() < tmp_event_time && tmp_event_time <= ti.end_height()) );

  // Sample the event type
  sampleEventType(tmp_event_time, tmp_event_line, ti, return_event);
}


/**
* Given that an event has happened, this function samples the events type. 
* 
* In particular requires that the 'temporary' forest members samples_, rates_, 
* active_nodes_, and nodes_timelines_ are set correctly beforehand.  
*/
void Forest::sampleEventType(const double &time, const size_t &time_line, 
                              const TimeInterval &ti, Event &event) const {
  event = Event(time);

  if ( rates_[time_line] == 0.0 ) throw std::logic_error("And event with rate 0 has happend!");

  // Situation where it is clear what happend:
  if (time == -1) return;
  if (time_line == 2) return event.setToCoalescence(active_node(1), 1);

  double sample = random_generator()->sample() * rates_[time_line];
  //std::cout << "Sample: " << sample << std::endl;

  for (size_t i = 0; i < 2; ++i) {
    // Only Nodes in state 1 or 2 can do something
    if ( states_[i] == 0 ) continue;

    // Coalescence can occur on all time lines  
    if (states_[i] == 1 && active_nodes_timelines_[i] == time_line) {
      sample -= calcCoalescenceRate(active_node(i)->population(), ti);
      //std::cout << i << ":" << "Sample after coalescence:" << sample << std::endl;
      if (sample <= 0.0) return event.setToCoalescence(active_node(i), i); 
    }

    // Migration and Recombination only on time line 0    
    if (time_line != 0) continue;

    // Recombination
    if (states_[i] == 2) {
      sample -= calcRecombinationRate(active_nodes_[i]);
      if (sample <= 0.0) return event.setToRecombination(active_node(i), i);
      continue;
    }

    // Migration
    assert( states_[i] == 1 );
    if ( sample < model().total_migration_rate(active_node(i)->population()) ) {
      //std::cout << i << ":" << "Migrating" << sample << std::endl;
      for ( size_t j = 0; j < model().population_number(); ++j) {
        sample -= model().migration_rate(active_node(i)->population(), j);
        if ( sample <= 0.0 ) return event.setToMigration(active_node(i), i, j);
      } 
      throw std::logic_error("Error Sampling Type of Event");
    }
    sample -= model().total_migration_rate(active_node(i)->population());
  }

  // If we are here, than we should have sampled a pw coalescence...
  if (! (states_[0] == 1 && states_[1] == 1 ))
    throw std::logic_error("Rates wrong!1");  
  if (! (active_nodes_[0]->population() == active_nodes_[1]->population() ))
    throw std::logic_error("Rates wrong!2");  
  if( sample > calcPwCoalescenceRate(active_nodes_[0]->population()) ) {
    throw std::logic_error("Rates wrong!3");  
  };
  return event.setToPwCoalescence();
}


/**
 * Looks if there was an event on a sampled time (e.g. this new_time != -1)
 * and if so sets current_time to the event with the smallest time and increases
 * the time_line counter.
 *
 * \param new_time An ExpoLimit Sample
 * \param time_line The timeline that the sample was from 
 * \param current_time The variable that save the time of the nearest event 
 * \param time_line The variable that saves the timeline of the nearest event 
 * \return Nothing, but updates current_time and current_time_line   
 */
void Forest::selectFirstTime(const double &new_time, const size_t &time_line, 
                             double &current_time, size_t &current_time_line) const {
  if (new_time == -1) return;
  if (current_time == -1 || new_time < current_time) {
    current_time = new_time;
    current_time_line = time_line;
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


double Forest::calcCoalescenceRate(const size_t &pop, const TimeInterval &ti) const {
    // Rate for each pair is 1/(2N), as N is the diploid population size
    return ( ti.numberOfContemporaries(pop) / ( 2.0 * this->model().population_size(pop) ) );
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
void Forest::implementCoalescence(const Event &event, TimeIntervalIterator &tii) {
  // Coalescence: sample target point and implement the coalescence
  assert( event.node() == active_node(event.active_node_nr()) );

  Node* coal_node = event.node();
  Node* target = (*tii).getRandomContemporary(coal_node->population());

  dout << "* * * Above node " << target << std::endl;
  assert( target->height() < event.time() ); 
  assert( coal_node->population() == target->population() );
  assert( getEventNode() != NULL );
  assert( getOtherNode() != NULL );

  Node* new_node;

  // Look if we can reuse the root that coalesced as new node
  if ( coal_node->numberOfChildren() == 1 && !coal_node->is_migrating() ){
    new_node = coal_node;
    coal_node = coal_node->first_child();
    nodes()->move(new_node, event.time());
  } else {
    new_node = new Node(event.time());
    new_node->change_child(NULL, coal_node);
    coal_node->set_parent(new_node);
    nodes()->add(new_node);
  }

  // Now we have:
  // target = Node in the target tree under the coalescence
  // coal_node =  Root of the coalescing tree 
  // new_node =   New parent of both other nodes

  // Update new_node
  new_node->set_population(coal_node->population());
  new_node->change_child(NULL, target);
  new_node->set_parent(target->parent());
  if (!target->local()) {
    new_node->set_local(false);
    new_node->set_last_update(target->last_update());
  } else {
    new_node->set_local(true);
  }

  // Integrate it into the tree
  target->set_parent(new_node);
  new_node->parent()->change_child(target, new_node);

  // And update
  updateAbove(coal_node, false, false); // Update coal_node ONLY
  // updateAbove(new_node, false, false); 
  // Updating the new_node will always activate it, but there may still be
  // unimplemented recombinations above

  set_active_node(event.active_node_nr(), new_node);

  if ( getOtherNodesState() == 2 ) {
    // If the coalescing node coalesced into the branch directly above 
    // a recombining node, then we are done.
    if ( getOtherNode()->parent() == getEventNode() ) {
      getOtherNode()->set_last_update(this->current_base());
      dout << "* * * Recombining Node moved into coalesced node. Done." << std::endl;
      updateAbove(getOtherNode(), false, false);
      updateAbove(getOtherNode()->parent());
      coalescence_finished_ = true;
      return;
    }

    // The branch below the event will be made local later anyway, so we don't
    // have to care about marking it as updated.
  }

  // Check if are can stop.
  if ( getEventNode()->local() ) {
    // Only one node should be active if there are still local nodes at the
    // current height. Hence we are done if so.
    assert( getOtherNodesState() == 0 );

    dout << "* * * We hit the local tree. Done." << std::endl;
    updateAbove(getEventNode()); 
    coalescence_finished_ = true;
    return;
  }

  if (active_node(0) == active_node(1)) {
    dout << "* * * Coalescend into other Active Node. Done." << std::endl;
    updateAbove(getEventNode()); 
    coalescence_finished_ = true;
    return;
  }

  // If we hit an non-local branch:
  // Begin next interval at the coalescence height and remove the branch
  // below from contemporaries.
  tii.splitCurrentInterval(getEventNode(), target);
}



/** 
 * @brief Modifies the forest to reflect that two coalescing nodes coalesced together.
 * 
 * @param root_1 The first coalescing node
 * @param root_2 The second coalescing node
 * @param time   The time at which the coalescence happens
 */
void Forest::implementPwCoalescence(Node* root_1, Node* root_2, const double &time) {
  dout << "* * Both nodes coalesced together" << std::endl;
  dout << "* * Implementing..." << std::flush;
  Node* new_root = NULL;
  assert( root_1 != NULL );
  assert( root_2 != NULL );
  assert( root_1->population() == root_2->population() );

  // both nodes may or may not mark the end of a single branch at the top of their tree,
  // which we don't need anymore.
  if (root_1->numberOfChildren() == 1 && !root_1->is_migrating()) {
    if (root_2->numberOfChildren() == 1 && !root_2->is_migrating()) {
      // both trees have a single branch => delete one
      root_2 = root_2->first_child();
      this->nodes()->remove(root_2->parent());
      root_2->set_parent(NULL);
      assert( root_2 != NULL );
      assert( root_2->is_root() );
    }
    // (now) only root_1 has a single branch => use as new root
    this->nodes()->move(root_1, time);
    new_root = root_1;
    root_1 = root_1->first_child();
    assert( root_1 != NULL );
  } 
  else if (root_2->numberOfChildren() == 1 && !root_2->is_migrating()) {
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
  new_root->set_population(root_1->population());

  updateAbove(root_1, false, false);
  updateAbove(root_2, false, false);
  updateAbove(new_root, false, false);
  dout << " done" << std::endl;

  assert( this->local_root()->height() == time );
}


void Forest::implementRecombination(const Event &event, TimeIntervalIterator &ti) {
  TreePoint event_point = TreePoint(event.node(), event.time(), false); //XXX beauty this up
  set_active_node(event.active_node_nr(), cut(event_point));
  updateAbove(event.node(), false, true); // Make the branch below the rec local
  //updateAbove(active_node(event.active_node_nr()), false, false);
  ti.recalculateInterval();

  assert( this->printTree() );
}


void Forest::implementMigration(const Event &event, TimeIntervalIterator &ti) {
  dout << "* * Implementing migration event... " << std::flush;
  assert( event.node()->is_root() );
  
  // There is only little to do if we can reuse the event.node() 
  if ( event.node()->is_unimportant() ) {
    dout << "Reusing: " << event.node() << "... " << std::flush;
    nodes()->move(event.node(), event.time());
    event.node()->set_population(event.mig_pop());
    updateAbove(event.node());
    assert(event.node()->is_migrating());
  }
  else {
    // Otherwise create a new node that marks the migration event
    Node* mig_node = new Node(event.time(), true);
    dout << "Marker: " << mig_node << "... " << std::flush; 
    nodes()->add(mig_node, event.node());
    mig_node->set_population(event.mig_pop());
    
    // Integrate it into the tree
    event.node()->set_parent(mig_node);
    mig_node->set_first_child(event.node());
    updateAbove(event.node(), false, false);
    updateAbove(mig_node);

    // And set it active
    this->set_active_node(event.active_node_nr(), mig_node);
    assert(mig_node->is_migrating());
  }
  // And recalculate the interval
  ti.recalculateInterval();
  dout << "done." << std::endl;
}


void Forest::implementFixedTimeEvent(TimeIntervalIterator &ti) {
  dout << "* * Fixed time event" << std::endl;
  for (size_t i = 0; i < 2; ++i) {
    if (states_[i] != 1) continue;
    double prob;
    for (size_t j = 0; j < model().population_number(); ++j) {
      prob = model().single_mig_pop(active_node(i)->population(), j);
      if (prob == 0.0) continue;
      if (prob == 1.0 || prob <= random_generator()->sample() ) {
        tmp_event_ = Event((*ti).start_height());
        tmp_event_.setToMigration(active_node(i), i, j);
        implementMigration(tmp_event_, ti);
      }
    }
  }
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
  assert ( node != NULL );
  if ( model().exact_window_length() == -1 ) return false;
  if ( node == local_root() ) return false;
  if ( node->is_migrating() ) return false;

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

bool areSame(const double &a, const double &b, const double &epsilon) {
  // from Knuths "The art of computer programming"
  return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


void Forest::initialize_recomb_coalescent(const double rec_height){};
void Forest::initialize_event(double start_time){};
void Forest::record_coalevent(const TimeInterval & current_event, double end_time){};
void Forest::record_recombevent(const TimeInterval & current_event, double end_time){};
void Forest::clear_initial_coalevent(){};

/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013, 2014 Paul R. Staab, Sha (Joe) Zhu, Dirk Metzler and Gerton Lunter
 * 
 * This file is part of scrm.
 * 
 * scrm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or * (at your option) any later version.  * * This program is distributed in the hope that it will be useful, * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "forest.h"

/******************************************************************
 * Debugging Utils
 *****************************************************************/

void Forest::createExampleTree() {
  this->clear();
  this->writable_model()->disable_approximation();
  // Only set the number of samples to 4, but keep rest of the model
  this->writable_model()->sample_times_.clear();
  this->writable_model()->sample_populations_.clear();
  this->writable_model()->addSampleSizes(0.0, std::vector<size_t>(1, 4));

  this->rec_bases_.push_back(5.0);  
  this->current_rec_ = 1;
  //this->rec_bases_.push_back(105.0);  

  Node* leaf1 = nodes()->createNode(0, 1);
  Node* leaf2 = nodes()->createNode(0, 2);
  Node* leaf3 = nodes()->createNode(0, 3);
  Node* leaf4 = nodes()->createNode(0, 4);

  leaf1->set_label(1);
  leaf2->set_label(2);
  leaf3->set_label(3);
  leaf4->set_label(4);

  this->nodes()->add(leaf4);
  this->nodes()->add(leaf3);
  this->nodes()->add(leaf2);
  this->nodes()->add(leaf1);

  Node* node12 = nodes()->createNode(1);
  this->addNodeToTree(node12, NULL, leaf1, leaf2);

  Node* node34 = nodes()->createNode(3);
  this->addNodeToTree(node34, NULL, leaf3, leaf4);

  Node* root = nodes()->createNode(10);
  this->addNodeToTree(root, NULL, node12, node34);
  this->set_local_root(root);
  this->set_primary_root(root);

  // Add a non-local tree
  Node* nl_node = nodes()->createNode(4); 
  nl_node->make_nonlocal(current_rec_);
  Node* nl_root = nodes()->createNode(6);
  nl_root->make_nonlocal(current_rec_);
  
  nl_node->set_parent(nl_root);
  nl_root->set_first_child(nl_node);
  this->nodes()->add(nl_node);
  this->nodes()->add(nl_root);
  updateAbove(nl_node);

  updateAbove(leaf1);
  updateAbove(leaf2);
  updateAbove(leaf3);
  updateAbove(leaf4);

  this->set_sample_size(4);

  this->contemporaries_ = ContemporariesContainer(model().population_number(), 
                                                  model().sample_size(),
                                                  random_generator());
  this->tmp_event_time_ = -1; 
  this->coalescence_finished_ = true;

  assert( this->checkTreeLength() );
  assert( this->checkTree() );
}

void Forest::createScaledExampleTree() {
  this->createExampleTree();

  this->nodes()->at(4)->set_height(1 * 4 * model().default_pop_size); 
  this->nodes()->at(5)->set_height(3 * 4 * model().default_pop_size); 
  this->nodes()->at(6)->set_height(4 * 4 * model().default_pop_size); 
  this->nodes()->at(7)->set_height(6 * 4 * model().default_pop_size); 
  this->nodes()->at(8)->set_height(10 * 4 * model().default_pop_size); 

  updateAbove(nodes()->at(4));
  updateAbove(nodes()->at(5));
  updateAbove(nodes()->at(6));

  assert( this->checkTreeLength() );
  assert( this->checkTree() );
}

double Forest::calcTreeLength() const {
  double local_length = 0;

  for (ConstNodeIterator it = getNodes()->iterator(); it.good(); ++it) {
    if ( *it == local_root() ) return local_length;
    if ( (*it)->is_root() || !(*it)->local() ) continue;
    local_length += (*it)->height_above();
  }

  return local_length;
}


void Forest::addNodeToTree(Node *node, Node *parent, Node *first_child, Node *second_child) {
  this->nodes()->add(node);

  if (parent != NULL) {
    node->set_parent(parent);
    if (parent->first_child() == NULL) parent->set_first_child(node);
    else {
      if (parent->first_child()->height() > node->height()) {
        parent->set_second_child(parent->first_child());
        parent->set_first_child(node);
      } else {
        parent->set_second_child(node);
      }
    }
  }

  if (first_child != NULL) {
    node->set_first_child(first_child);
    first_child->set_parent(node);
  }

  if (second_child != NULL) {
    node->set_second_child(second_child);
    second_child->set_parent(node);
  }
}


bool Forest::checkTreeLength() const {
  double local_length = calcTreeLength();

  if ( !areSame(local_length, getLocalTreeLength(), 0.000001) ) {
    scrmdout << "Error: local tree length is " << this->getLocalTreeLength() << " ";
    dout << "but should be " << local_length << std::endl;
    return(0);
  }

  return(1);
}


bool Forest::checkInvariants(Node const* node) const {
  if (node == NULL) {
    bool okay = 1;

    for (ConstNodeIterator it = getNodes()->iterator(); it.good(); ++it) {
      if ( (*it)->height() >= local_root()->height()) {
        if (!(*it)->local()) continue;
        dout << "Node " << *it << " is above the local root and local!" << std::endl;
        okay = 0;  
      } else {
        okay *= checkInvariants(*it);
      }
    }
    return(okay);
  }

  size_t samples_below = node->in_sample();
  double length_below = 0;

  if (node->first_child() != NULL) {
    samples_below += node->first_child()->samples_below();
    length_below += node->first_child()->length_below();
    if (node->first_child()->local()) 
      length_below += node->first_child()->height_above();
  }

  if (node->second_child() != NULL) {
    samples_below += node->second_child()->samples_below();
    length_below += node->second_child()->length_below();
    if (node->second_child()->local()) 
      length_below += node->second_child()->height_above();
  }

  if ( samples_below != node->samples_below() ||
      !areSame(length_below, node->length_below(), 0.00001) ) {
    dout << "Node " << node << " not up to date" << std::endl;
    dout << "samples_below: is " << node->samples_below() 
         << " and should be " << samples_below << std::endl;
    dout << "length_below: is " << node->length_below() 
         << " and should be " << length_below 
         << " ( Diff " << node->length_below() - length_below << " )" << std::endl;

    printNodes();
    printTree();
    return false;
  }

  if ( (samples_below == 0 || samples_below == sample_size()) && node->local() ) {
    dout << "Node " << node << " is local but should be non-local" << std::endl;
    return false;
  }

  return true;
}


bool Forest::checkLeafsOnLocalTree(Node const* node) const {
  if (node == NULL) {
    size_t all_on_tree = 1;
    bool on_tree = 0;
    for (ConstNodeIterator it = getNodes()->iterator(); it.good(); ++it) {
      if ( !(*it)->in_sample() ) continue;
      on_tree = checkLeafsOnLocalTree(*it);
      if (!on_tree) dout << "Leaf " << *it << " is not on local tree!" << std::endl;
      all_on_tree *= on_tree;
    }
    return(all_on_tree);
  }
  if ( node->local() ) return( checkLeafsOnLocalTree(node->parent()) );
  return( node == this->local_root() );
}


bool Forest::checkNodeProperties() const {
  bool success = true;
  for (ConstNodeIterator it = getNodes()->iterator(); it.good(); ++it) {
    if ( !(*it)->local() ) {
      if ( (*it)->last_update() == 0 && !(*it)->is_root() ) {
        dout << "Error: Node " << *it << " non-local without update info" << std::endl;
        success = false;
      }
    }
  } 
  return success; 
} 


bool Forest::checkTree(Node const* root) const {
  if (root == NULL) {
    bool good = true;
    // Default when called without argument
    for (ConstNodeIterator it = getNodes()->iterator(); it.good(); ++it) {
      if ( (*it)->is_root() ) good *= checkTree(*it);
    }

    good *= this->checkInvariants();
    good *= this->checkNodeProperties();
    good *= this->checkTreeLength();
    good *= this->checkRoots();
    return good;
  }
  assert( root != NULL );

  Node* h_child = root->second_child();
  Node* l_child = root->first_child();

  bool child1 = 1;
  if (h_child != NULL) {
    if (l_child == NULL) {
      dout << root << ": only child is second child" << std::endl;
      return 0;
    }
    if (h_child->parent() != root) {
      dout << h_child << ": is child of non-parent" << std::endl;
      return 0;
    }
    if (h_child->height() > root->height()) {
      dout << root << ": has child with greater height" << std::endl;
      return 0;
    }
    if (h_child->population() != root->population()) {
      dout << root << ": has child of other population" << std::endl;
      return 0;
    }
    if (l_child->population() != root->population()) {
      dout << root << ": has child of other population" << std::endl;
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

    if (l_child->height() > root->height()) {
      dout << root << ": has child with greater height" << std::endl;
      return 0;
    }
  }

  // Check that parent if above node
  if (!root->is_root()) {
    Node* parent = root->parent();
    Node const* current = root;
    while (current != parent) {
      if (current->is_last()) {
        dout << root << ": node is above it's parent.";
        return 0;
      }
      current = current->next();
    } 
  }

  return child1*child2;
}




/******************************************************************
 * Tree Printing
 *****************************************************************/
bool Forest::printTree() const {
  //this->printNodes();
  std::vector<Node const*> positions = this->determinePositions();
  //this->printPositions(positions);
  std::vector<Node const*>::iterator position;
  int h_line;
  double start_height = 0, 
         end_height = getNodes()->get(0)->height();

  for (ConstNodeIterator ni = getNodes()->iterator(); ni.good(); ) {
    if ( !(*ni)->is_root() && (*ni)->height_above() == 0.0 ) {
      std::cout << "A rare situation occurred were a parent and a child have exactly "
           << "the same height. We can't print such trees here, the algorithm however"
           << "should not be affected." << std::endl; 
      return 1;
    }
    h_line = 0;
    start_height = end_height;
    while ( ni.height() <= end_height ) ++ni;
    end_height = ni.height(); 
    //std::cout << start_height << " - " << end_height << std::endl;

    for (position = positions.begin(); position != positions.end(); ++position) {
      assert( *position != NULL );
      if ( (*position)->height() == start_height ) {
        if ( (*position)->local() || *position == local_root() ) dout << "╦";
        else dout << "┬";
        if ( (*position)->numberOfChildren() == 2 ) {
          h_line = 1 + !((*position)->local());
          if ( *position == local_root() ) h_line = 1;
        }
        if ( (*position)->numberOfChildren() == 1 ) {
          h_line = 0;
        }
      } 
      else if ( (*position)->height() < start_height &&
                (*position)->parent_height() >= end_height ) {
        if ( (*position)->local() ) dout << "║";
        else dout << "│";

      } 
      else if ( (*position)->parent_height() == start_height ) {
        if ( *position == (*position)->parent()->first_child() ) {
          if ( (*position)->local() ) { 
            dout << "╚";
            h_line = 1;
          }
          else {
            dout << "└";
            h_line = 2;
          }
        }
        else {
          if ( (*position)->local() ) dout << "╝";
          else dout << "┘";
          h_line = 0;
        }
      }
      else {
        if ( h_line == 0 ) dout << " ";
        else if ( h_line == 1 ) dout << "═";
        else dout << "─";
      }
    }
    dout << " - " << std::setw(7) << setprecision(7) << std::right << start_height << " - "; 
    for (position = positions.begin(); position != positions.end(); ++position) {
      if (*position == NULL) continue;
      if ( (*position)->height() == start_height ) {
        if ((*position)->label() != 0) dout << (*position)->label() << ":";
        if (!(*position)->is_migrating()) dout << *position << "(" << (*position)->population() << ") ";
        else dout << *position << "(" << (*position)->first_child()->population()
                  << "->" << (*position)->population() << ") ";
        if (nodeIsOld(*position)) dout << "old ";
      }
    }
    dout << std::endl;
  }
  return true;
}


void Forest::printTree_cout() const {
  //this->printNodes();
  std::vector<Node const*> positions = this->determinePositions();
  //this->printPositions(positions);
  std::vector<Node const*>::iterator position;
  int h_line;
  double start_height = 0, 
         end_height = getNodes()->get(0)->height();

  for (ConstNodeIterator ni = getNodes()->iterator(); ni.good(); ) {
    if ( !(*ni)->is_root() && (*ni)->height_above() == 0.0 ) {
      std::cout << "A rare situation occurred were a parent and a child have exactly "
           << "the same height. We can't print such trees here, the algorithm however"
           << "should not be affected." << std::endl; 
      //return 1;
    }
    h_line = 0;
    start_height = end_height;
    while ( ni.height() <= end_height ) ++ni;
    end_height = ni.height(); 
    //std::cout << start_height << " - " << end_height << std::endl;

    for (position = positions.begin(); position != positions.end(); ++position) {
      assert( *position != NULL );
      if ( (*position)->height() == start_height ) {
        if ( (*position)->local() || *position == local_root() ) std::cout << "╦";
        else std::cout << "┬";
        if ( (*position)->countChildren() == 2 ) {
          h_line = 1 + !((*position)->local());
          if ( *position == local_root() ) h_line = 1;
        }
        if ( (*position)->countChildren() == 1 ) {
          h_line = 0;
        }
      } 
      else if ( (*position)->height() < start_height &&
                (*position)->parent_height() >= end_height ) {
        if ( (*position)->local() ) std::cout << "║";
        else std::cout << "│";

      } 
      else if ( (*position)->parent_height() == start_height ) {
        if ( *position == (*position)->parent()->first_child() ) {
          if ( (*position)->local() ) { 
            std::cout << "╚";
            h_line = 1;
          }
          else {
            std::cout << "└";
            h_line = 2;
          }
        }
        else {
          if ( (*position)->local() ) std::cout << "╝";
          else std::cout << "┘";
          h_line = 0;
        }
      }
      else {
        if ( h_line == 0 ) std::cout << " ";
        else if ( h_line == 1 ) std::cout << "═";
        else std::cout << "─";
      }
    }
    std::cout << " - " << std::setw(7) << std::setprecision(7) << std::right << start_height << " - "; 
    for (position = positions.begin(); position != positions.end(); ++position) {
      if (*position == NULL) continue;
      if ( (*position)->height() == start_height ) {
        if ((*position)->label() != 0) std::cout << (*position)->label() << ":";
        if (!(*position)->is_migrating()) std::cout << *position << "(" << (*position)->population() << ") ";
        else std::cout << *position << "(" << (*position)->first_child()->population()
                  << "->" << (*position)->population() << ") ";
        if (nodeIsOld(*position)) std::cout << "old ";
      }
    }
    std::cout << std::endl;
  }
}


/**
 *  For printing the tree, each node gets assigned its own column in the printed area, 
 *  referred to as its positions. This function determines the position for all
 *  nodes and returns the nodes in a vector sorted by position.   
 *
 *  \return Vector of all nodes, sorted by position 
 */
std::vector<Node const*> Forest::determinePositions() const {
  std::vector<Node const*> positions(this->getNodes()->size(), NULL); 

  ReverseConstNodeIterator it;
  std::vector<const Node*>::iterator cit;
  size_t lines_left, lines_right, position, root_offset = 0;
  Node const* current_node;

  for (it = getNodes()->reverse_iterator(); it.good(); ++it) {
    current_node = *it;

    lines_left = countLinesLeft(current_node);
    lines_right = countLinesRight(current_node);

    if ( current_node->is_root() ) {
      // Add root to the right of all current trees
      position = countBelowLinesLeft(current_node->first_child()) + lines_left + root_offset;
      //std::cout << current_node << " " << position << " " << lines_left << " "
      //          << lines_right << " "  
      //          << countBelowLinesLeft(current_node->first_child()) << std::endl;
      
      root_offset = position + 
                    countBelowLinesRight(current_node->second_child()) + 
                    lines_right + 1;

      assert( positions[position] == NULL ); 
      positions[position] = current_node;
    } else {
      // Get the position of the node (which was assigned when looking at its
      // parent
      position = 0;
      for (cit = positions.begin(); cit < positions.end(); ++cit) {
        if ( *cit == current_node ) break;
        ++position;
      }
    }
    
    // Insert the child/children into branches
    if (current_node->first_child() != NULL) {
        assert( positions.at(position - lines_left) == NULL );         
        positions[position - lines_left] =  current_node->first_child();
    }
    

    if (current_node->second_child() != NULL) { 
        assert( positions.at(position + lines_right) == NULL );         
        positions[position + lines_right] = current_node->second_child();
    }
  
  }
  return positions;
}

  void Forest::printPositions(const std::vector<Node const*> &positions) const {
      for (size_t col = 0; col < positions.size() ; ++col) {
        std::cout << positions[col] << " ";
      } 
      std::cout << std::endl;
  }

  int Forest::countLinesLeft(Node const* node) const {
    if ( node->first_child() == NULL ) return 0;
    //if ( node->second_child() == NULL ) return 1;
    return ( 1 + countBelowLinesRight(node->first_child()) );
  }

  int Forest::countLinesRight(Node const* node) const {
    if ( node->first_child() == NULL ) return 0;
    if ( node->second_child() == NULL ) return 0;
    return ( 1 + countBelowLinesLeft(node->second_child()) );
  }

  int Forest::countBelowLinesLeft(Node const* node) const {
    if ( node == NULL ) return 0;
    if ( node->first_child() == NULL ) return 0;
    else return ( countLinesLeft(node) + countBelowLinesLeft(node->first_child()) );
  }

  int Forest::countBelowLinesRight(Node const* node) const {
    if ( node == NULL ) return 0;
    if ( node->second_child() == NULL ) return 0;
    else return ( countLinesRight(node) + countBelowLinesRight(node->second_child()) );
  }

  bool Forest::printNodes() const {
    scrmdout << std::setw(10) << std::right << "Node";
    dout << std::setw(10) << std::right << "Height";
    dout << std::setw(6) << std::right << "label";
    dout << std::setw(10) << std::right << "Parent";
    dout << std::setw(10) << std::right << "1th_child";
    dout << std::setw(10) << std::right << "2nd_child";
    dout << std::setw(6) << std::right << "local";
    dout << std::setw(6) << std::right << "pop";
    dout << std::setw(10) << std::right << "l_upd";
    dout << std::setw(6) << std::right << "s_bel";
    dout << std::setw(10) << std::right << "l_bel";
    dout << std::endl;

    for(size_t i = 0; i < this->getNodes()->size(); ++i) {
      scrmdout << std::setw(10) << std::right << this->getNodes()->get(i);
      dout << std::setw(10) << std::right << this->getNodes()->get(i)->height();
      dout << std::setw(6) << std::right << this->getNodes()->get(i)->label();
      if (!getNodes()->get(i)->is_root()) 
        dout << std::setw(15) << std::right << this->getNodes()->get(i)->parent();
      else dout << std::setw(15) << std::right << 0;
      dout << std::setw(15) << std::right << this->getNodes()->get(i)->first_child();
      dout << std::setw(15) << std::right << this->getNodes()->get(i)->second_child();
      dout << std::setw(6) << std::right << this->getNodes()->get(i)->local();
      dout << std::setw(6) << std::right << this->getNodes()->get(i)->population();
      dout << std::setw(10) << std::right << this->getNodes()->get(i)->last_update();
      dout << std::setw(6) << std::right << this->getNodes()->get(i)->samples_below();
      dout << std::setw(10) << std::right << this->getNodes()->get(i)->length_below();
      dout << std::endl;
    }
    scrmdout << "Local Root:    " << this->local_root() << std::endl;
    scrmdout << "Primary Root:  " << this->primary_root() << std::endl;
    return true;
  }


bool Forest::checkForNodeAtHeight(const double height) const {
  for (auto it = getNodes()->iterator(); it.good(); ++it) {
    if ((*it)->height() == height) return true;
    if ((*it)->height() > height) return false;
  }
  return false;
}

// Checks if all nodes in contemporaries are contemporaries.
bool Forest::checkContemporaries(const double time) const {
  // Check if all nodes in contemporaries() are contemporaries
  for (size_t pop = 0; pop < model().population_number(); ++pop) {
    for (auto it = contemporaries_.begin(pop); it != contemporaries_.end(pop); ++it) {
      if ( *it == NULL ) {
        dout << "NULL in contemporaries" << std::endl; 
        return 0;
      }

      if ( (*it)->is_root() ) {
        dout << "Root " << *it << " in contemporaries" << std::endl; 
        return 0;
      }

      if ( (*it)->height() > time || (*it)->parent_height() <= time ) {
        dout << "Non-contemporary node " << *it << " in contemporaries " 
             << "at time " << time << " (node at " << (*it)->height() 
             << "; parent at " << (*it)->parent_height() << ")." << std::endl; 
        printNodes();
        return 0;
      }

      if ( nodeIsOld(*it) ) { 
        if ( *it == local_root() ) {
          if ( !(*it)->is_root() ) {
            dout << "Branch above local root should be pruned but is not" << std::endl;
            return 0;
          }
        } else {
          dout << "Contemporary node " << *it << " should be pruned by now!" << std::endl;
          return 0;
        }
      }

      for (size_t i = 0; i < 2; ++i) {
        if ( *it == active_node(i) && states_[i] == 1 ) {
          dout << "Coalescing node a" << i << " in contemporaries!" << std::endl;
          return 0;
        }
      }
    }
  }

  // Check if all contemporaries are in contemporaries() 
  for (auto ni = getNodes()->iterator(); ni.good(); ++ni) {
    if ( (*ni)->height() <= time && time < (*ni)->parent_height()) {
      if ( *ni == active_node(0) && states_[0] == 1 ) continue; 
      if ( *ni == active_node(1) && states_[1] == 1 ) continue; 
      
      bool found = false;
      size_t pop = (*ni)->population();
      for (auto it = contemporaries_.begin(pop); it != contemporaries_.end(pop); ++it) {
        if ( *it == *ni ) {
          found = true;
          break;
        }
      } 
      if (!found) { 
        scrmdout << "Node " << *ni << " (height " << (*ni)->height() 
             <<") not in contemporaries." << std::endl;
        dout << ti.start_height() << " " << ti.end_height() << std::endl;
        return 0;
      }
    }
  }

  return 1;
}

bool Forest::checkRoots() const {
  // Check that local_root() really is the local root:
  if (local_root()->samples_below() != sample_size() ||
      local_root()->first_child() == NULL ||
      local_root()->second_child() == NULL ||
      (!local_root()->first_child()->local()) ||
      (!local_root()->second_child()->local()) ) {
    dout << local_root() << " is registered as local root, but is not." << std::endl;
    return false;
  } 
  
  // Check that primary_root() really is the primary root:
  Node* node = local_root();
  while (!node->is_root()) node = node->parent(); 
  if (node != primary_root()) {
    dout << primary_root() << " is registered as primary root, but "
         << node << " is." << std::endl;
    return false;
  }

  return true;
}

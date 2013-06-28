#include "forest.h"
#include <sstream>

/******************************************************************
 * Debugging Utils
 *****************************************************************/

void Forest::createExampleTree() {
  this->nodes()->clear();
  Node* leaf1 = new Node(0, true, 0, 1);
  Node* leaf2 = new Node(0, true, 0, 1);
  Node* leaf3 = new Node(0, true, 0, 1);
  Node* leaf4 = new Node(0, true, 0, 1);
  this->nodes()->add(leaf1);
  this->nodes()->add(leaf2);
  this->nodes()->add(leaf3);
  this->nodes()->add(leaf4);


  Node* node12 = new Node(1);
  this->addNodeToTree(node12, NULL, leaf1, leaf2);

  Node* node34 = new Node(3);
  this->addNodeToTree(node34, NULL, leaf3, leaf4);

  Node* root = new Node(10);
  this->addNodeToTree(root, NULL, node12, node34);
  this->set_local_root(root);
  this->set_primary_root(root);

  // Add a non-local tree
  Node* nl_node = new Node(4, false, 5, 0); 
  Node* nl_root = new Node(6, false, 5, 0);
  nl_node->set_parent(nl_root);
  nl_root->set_first_child(nl_node);
  this->nodes()->add(nl_node);
  this->nodes()->add(nl_root);
  updateAbove(nl_node);

  updateAbove(node12);
  updateAbove(node34);

  this->set_current_base(5);  
  assert( this->checkTreeLength() );
  assert( this->checkTree() );
}


double Forest::calcTreeLength() const {
  double local_length = 0;

  for (ConstNodeIterator it = getNodes()->iterator(); it.good(); ++it) {
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

  if ( !areSame(local_length, local_tree_length()) ) {
    dout << "Error: local tree length is " << this->local_tree_length() << " ";
    dout << "but should be " << local_length << std::endl;
    return(0);
  }

  return(1);
}


bool Forest::checkInvariants(Node const* node) const {
  if (node == NULL) {
    bool okay = 1;

    for (ConstNodeIterator it = getNodes()->iterator(); it.good(); ++it) {
      okay *= checkInvariants(*it);
    }
    return(okay);
  }

  size_t samples_below = node->in_sample();
  double length_below = 0;

  // Allow non-local roots of have false invariants
  //if ( node->is_root() && !node->local() ) return true;

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
  if ( !node->is_root() ) return( checkLeafsOnLocalTree(node->parent()) );
  return( node == this->primary_root() );
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
    //Default when called without argument
    for (ConstNodeIterator it = getNodes()->iterator(); it.good(); ++it) {
      if ( (*it)->is_root() ) good = checkTree(*it);
    }

    assert( this->checkInvariants() );
    assert( this->checkNodeProperties() );
    assert( this->checkLeafsOnLocalTree() );
    return good;
  }
  assert( root != NULL);

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




/******************************************************************
 * Tree Printing
 *****************************************************************/
bool Forest::printTree() {
  //this->printNodes();
  std::vector<Node const*> positions = this->determinePositions();
  //this->printPositions(positions);
  std::vector<Node const*>::iterator position;
  int h_line;

  for (TimeIntervalIterator tii = TimeIntervalIterator(this, getNodes()->get(0), false); 
       tii.good(); ++tii) {
    h_line = 0;
    for (position = positions.begin(); position != positions.end(); ++position) {
      assert( *position != NULL );
      if ( (*position)->height() == (*tii).start_height() ) {
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
      else if ( (*position)->height() < (*tii).start_height() &&
                (*position)->parent_height() >= (*tii).end_height() ) {
        if ( (*position)->local() ) dout << "║";
        else dout << "│";

      } 
      else if ( (*position)->parent_height() == (*tii).start_height() ) {
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
    dout << " - " << std::setw(7) << setprecision(7) << std::right << (*tii).start_height() << " - "; 
    for (position = positions.begin(); position != positions.end(); ++position) {
      if (*position == NULL) continue;
      if ( (*position)->height() == (*tii).start_height() ) {
        dout << *position << " ";
      }
    }
    dout << std::endl;
  }
  return true;
}

/******************************************************************
 * Tree Printing
 *****************************************************************/
bool Forest::printTree_cout() {
  //this->printNodes();
  std::vector<Node const*> positions = this->determinePositions();
  //this->printPositions(positions);
  std::vector<Node const*>::iterator position;
  int h_line;

  for (TimeIntervalIterator tii = TimeIntervalIterator(this, getNodes()->get(0), false); 
       tii.good(); ++tii) {
    h_line = 0;
    for (position = positions.begin(); position != positions.end(); ++position) {
      assert( *position != NULL );
      if ( (*position)->height() == (*tii).start_height() ) {
        if ( (*position)->local() || *position == local_root() ) cout << "╦";
        else cout << "┬";
        if ( (*position)->numberOfChildren() == 2 ) {
          h_line = 1 + !((*position)->local());
          if ( *position == local_root() ) h_line = 1;
        }
        if ( (*position)->numberOfChildren() == 1 ) {
          h_line = 0;
        }
      } 
      else if ( (*position)->height() < (*tii).start_height() &&
                (*position)->parent_height() >= (*tii).end_height() ) {
        if ( (*position)->local() ) cout << "║";
        else cout << "│";

      } 
      else if ( (*position)->parent_height() == (*tii).start_height() ) {
        if ( *position == (*position)->parent()->first_child() ) {
          if ( (*position)->local() ) { 
            cout << "╚";
            h_line = 1;
          }
          else {
            cout << "└";
            h_line = 2;
          }
        }
        else {
          if ( (*position)->local() ) cout << "╝";
          else cout << "┘";
          h_line = 0;
        }
      }
      else {
        if ( h_line == 0 ) cout << " ";
        else if ( h_line == 1 ) cout << "═";
        else cout << "─";
      }
    }
    cout << " - " << std::setw(7) << setprecision(7) << std::right << (*tii).start_height() << " - "; 
    for (position = positions.begin(); position != positions.end(); ++position) {
      if (*position == NULL) continue;
      if ( (*position)->height() == (*tii).start_height() ) {
        cout << *position << " ";
      }
    }
    cout << std::endl;
  }
  return true;
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
    dout << std::setw(10) << std::right << "Node";
    dout << std::setw(10) << std::right << "Height";
    dout << std::setw(10) << std::right << "Parent";
    dout << std::setw(10) << std::right << "h_child";
    dout << std::setw(10) << std::right << "l_child";
    dout << std::setw(6) << std::right << "local";
    dout << std::setw(6) << std::right << "l_upd";
    dout << std::setw(6) << std::right << "s_bel";
    dout << std::setw(10) << std::right << "l_bel";
    dout << std::endl;

    //   \t\t Height\t Parent\t h_child\t  l_child" << std::endl;
    //    std::cout << std::setw(20) << std::right << "Hi there!" << std::endl;
    //        std::cout << std::setw(20) << std::right << "shorter" << std::endl;
    for(size_t i = 0; i < this->getNodes()->size(); ++i) {
      dout << std::setw(10) << std::right << this->getNodes()->get(i);
      dout << std::setw(10) << std::right << this->getNodes()->get(i)->height();
      dout << std::setw(10) << std::right << this->getNodes()->get(i)->parent();
      dout << std::setw(10) << std::right << this->getNodes()->get(i)->second_child();
      dout << std::setw(10) << std::right << this->getNodes()->get(i)->first_child();
      dout << std::setw(6) << std::right << this->getNodes()->get(i)->local();
      dout << std::setw(6) << std::right << this->getNodes()->get(i)->last_update();
      dout << std::setw(6) << std::right << this->getNodes()->get(i)->samples_below();
      dout << std::setw(10) << std::right << this->getNodes()->get(i)->length_below();
      if ( this->getNodes()->get(i) != getNodes()->last() )
        dout << std::setw(10) << std::right << this->getNodes()->get(i)->next();
      dout << std::endl;
    }
    dout << "Local Root:    " << this->local_root() << std::endl;
    dout << "Primary Root:  " << this->primary_root() << std::endl;
    return true;
  }

  bool areSame(const double &a, const double &b, const double &epsilon) {
    // from Knuths "The art of computer programming"
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
  }

#include "forest.h"

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
  nl_root->set_lower_child(nl_node);
  this->nodes()->add(nl_node);
  this->nodes()->add(nl_root);
  updateAbove(nl_node);

  updateAbove(node12);
  updateAbove(node34);
  
  assert( this->checkTreeLength() );
  assert( this->checkTree() );
}


double Forest::calcTreeLength() const {
  double local_length = 0;

  for (ConstNodeIterator it = nodes()->iterator(); it.good(); ++it) {
    if ( (*it)->is_root() || !(*it)->local() ) continue;
    local_length += (*it)->height_above();
  }
  return local_length;
}


void Forest::addNodeToTree(Node *node, Node *parent, Node *lower_child, Node *higher_child) {
  this->nodes()->add(node);

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
  }

  if (lower_child != NULL) {
    node->set_lower_child(lower_child);
    lower_child->set_parent(node);
  }

  if (higher_child != NULL) {
    node->set_higher_child(higher_child);
    higher_child->set_parent(node);
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


bool Forest::checkLeafsOnLocalTree(Node const* node) const {
  if (node == NULL) {
    size_t all_on_tree = 1;
    bool on_tree = 0;
    for (ConstNodeIterator it = nodes()->iterator(); it.good(); ++it) {
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
  for (ConstNodeIterator it = nodes()->iterator(); it.good(); ++it) {
    if ( !(*it)->local() ) {
      if ( (*it)->last_update() == 0 ) {
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
    for (ConstNodeIterator it = nodes()->iterator(); it.good(); ++it) {
      if ( (*it)->is_root() ) good = checkTree(*it);
    }

    assert( this->checkNodeProperties() );
    return good;
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
    //This is only allowed for root and fake nodes
    if ( !root->is_root() ) { 
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




/******************************************************************
 * Tree Printing
 *****************************************************************/
bool Forest::printTree() const {
  ReverseConstNodeIterator it;
  std::vector<const Node*>::iterator cit;
  Node const* current_node;
  size_t lines_left, lines_right, position, root_offset;
  //double current_height = FLT_MAX;
  double current_height = -1;

  std::vector<const Node*> branches;
  for (size_t i=0; i < this->getNodes()->size(); ++i) branches.push_back(NULL);
  root_offset = 0;

  for (it = getNodes()->reverse_iterator(); it.good(); ++it) {
    current_node = *it;

    lines_left = countLinesLeft(current_node);
    lines_right = countLinesRight(current_node);

    if ( current_node->is_root() ) {
      // Add root to the right of all current trees
      position = countBelowLinesLeft(current_node) + lines_left + root_offset;
      root_offset = position + countBelowLinesRight(current_node) + lines_right + 1;
      branches[position] = current_node;
    } else {
      // Get the position of the node (which was assigned when looking at is
      // parent
      position = 0;
      for (cit = branches.begin(); cit < branches.end(); ++cit) {
        if ( *cit == current_node ) break;
        ++position;
      }
    }

    // Insert the child/children into branches
    if (current_node->lower_child() != NULL) {
      if (current_node->higher_child() == NULL){ 
        branches[position + 1] = current_node->lower_child();
      } else {
        branches[position - lines_left] =  current_node->lower_child();
        branches[position + lines_right] = current_node->higher_child();
      }
    }
   
    // Print branches
    /* 
    for (std::vector<const Node*>::iterator bit = branches.begin(); 
         bit != branches.end(); ++bit) {
        if (*bit == NULL) std::cout << "0        ";
        else std::cout << *bit << " ";
    }
    std::cout << std::endl;
    continue;
    */


    if (areSame(current_height, current_node->height())) {
      dout << current_node << " ";
      continue;
    }
    dout << std::endl;

    current_height = current_node->height();
    std::string nodes = "";

    for (size_t i=0; i < branches.size(); ++i) {
      
      if ( i < position - lines_left) {
        if ( branches[i] == NULL ) dout << " ";
        else if (areSame(branches[i]->height(), current_node->height())) 
          if (branches[i]->local()) dout << "+";
          else dout << "°";
        else if ( branches[i]->is_root() ) dout << " ";
        else dout << "|";
      } 
      else if ( i < position) dout << "-";
      else if ( i == position) { 
        if (current_node->local()) dout << "+";
        else dout << "°";
      }
      else if ( i <= position + lines_right) dout << "-";
      else {
        if ( branches[i] == NULL ) dout << " ";
        else if (areSame(branches[i]->height(), current_node->height())) 
          if (branches[i]->local()) dout << "+";
          else dout << "°";
        else if ( branches[i]->is_root() ) dout << " ";
        else dout << "|";
      }
      
      // remove Nodes after they were printed
      if ( branches[i] != NULL && areSame(branches[i]->height(), current_height) ) {
        branches[i] = NULL;
      }
    }
    dout << " ---" << current_node->height() << "\t" << current_node << " "; 
  }
  dout << std::endl;
  dout << "(Nodes on the right are not sorted)" << std::endl << std::endl;
  //Returns true, so that it can be place in asserts
  return true;
}

int Forest::countLinesLeft(Node const* node) const {
  if ( node->lower_child() == NULL ) return 0;
  if ( node->higher_child() == NULL ) return 0;
  return ( 1 + countBelowLinesRight(node) );
}

int Forest::countLinesRight(Node const* node) const {
  if ( node->lower_child() == NULL ) {
    if ( node->higher_child() == NULL ) return 0;
    return 1;
  }
  return ( 1 + countBelowLinesLeft(node) );
}

int Forest::countBelowLinesLeft(Node const* node) const {
  if ( node->lower_child() == NULL ) return 0;
  else return ( countLinesLeft(node->lower_child()) + countBelowLinesLeft(node->lower_child()) );
}

int Forest::countBelowLinesRight(Node const* node) const {
  if ( node->higher_child() == NULL ) return 0;
  else return ( countLinesRight(node->higher_child()) + countBelowLinesRight(node->higher_child()) );
}

bool Forest::printNodes() const {
  for(size_t i = 0; i < this->getNodes()->size(); ++i) {
    dout << "Addr: " << this->getNodes()->get(i) << " | ";
    dout << "Height: " << this->getNodes()->get(i)->height() << " | ";
    dout << "Parent: " << this->getNodes()->get(i)->parent() << " | ";
    dout << "h_child: " << this->getNodes()->get(i)->higher_child() << " | ";
    dout << "l_child: " << this->getNodes()->get(i)->lower_child() << std::endl;
  }
  dout << "Local Root:    " << this->local_root() << std::endl;
  dout << "Primary Root:  " << this->primary_root() << std::endl;
  return true;
}

bool areSame(double a, double b) {
  return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
  //return std::fabs(a - b) < 0.0001;
}

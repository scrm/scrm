#include "node_container.h"

/*******************************************************
 * Con- & Destructor
 *******************************************************/

NodeContainer::NodeContainer() {
  set_first(NULL);
  set_last(NULL);
  unsorted_node_ = NULL;
};


/*******************************************************
 * Iterators
 *******************************************************/

NodeIterator NodeContainer::iterator() { 
  return(NodeIterator(*this)); 
}


ConstNodeIterator NodeContainer::iterator() const { 
  return(ConstNodeIterator(*this)); 
}


ReverseConstNodeIterator NodeContainer::reverse_iterator() const { 
  return(ReverseConstNodeIterator(*this)); 
}



/*******************************************************
 * Management of Nodes
 *******************************************************/

Node* NodeContainer::at(size_t nr) const {
  Node* current = first();

  for (size_t i=0; i < nr; ++i) {
    assert(current != NULL);
    current = current->next();
  }

  if ( current == NULL ) throw std::out_of_range("NodeContainer out of range"); 
  return current;
}


// Adds 'node' to the container
// If you know that the node is higher in the tree than a other node,
// than you can specify the latter as 'after_node' to speedup the process.
void NodeContainer::add(Node* node, Node* after_node) {
  if (first() == NULL) {
    this->set_first(node);
    this->set_last(node);
    return;
  }
  assert(first() != NULL);

  // Before first node?
  if (node->height() <= first()->height()) {
    node->set_next(first());
    node->set_previous(NULL);
    first()->set_previous(node);
    this->set_first(node);
    assert( this->sorted() );
    return;
  }

  // After last node?
  if ( node->height() >= last()->height() ) {
    node->set_previous(last());
    node->set_next(NULL);
    last()->set_next(node);
    this->set_last(node);
    assert( this->sorted() );
    return;
  }

  assert( after_node == NULL || node->height() >= after_node->height() );

  if (after_node == NULL) after_node = first();
  Node* current = after_node;
  // Find position in between
  while ( current->height() <= node->height() ) {
    assert( !current->is_last() );
    if ( (!current->is_root()) && current->parent_height() < node->height() ) 
      current = current->parent();
    else current = current->next();
  }
 
  // And add the node;
  this->add_before(node, current);
  assert( this->sorted() );
};


void NodeContainer::remove(Node *node, const bool &del) {
  if ( node->is_first() && node->is_last() ) {
    this->set_first(NULL);
    this->set_last(NULL);
    return;
  }

  if ( node->is_first() ) {
    this->set_first(node->next());
    node->next()->set_previous(NULL);
    return;
  }

  if ( node->is_last() ) {
    this->set_last(node->previous());
    node->previous()->set_next(NULL);
    return;
  }

  node->previous()->set_next(node->next());
  node->next()->set_previous(node->previous());
  if (del){ delete node;}
  assert( this->sorted() );
}


void NodeContainer::move(Node *node, const double &new_height) {
  assert( node != NULL );

  // Stupid edge case first: We may have only one node.
  if ( node->is_first() && node->is_last() ) {
    node->set_height(new_height);
    return;
  }

  // Remove from old place
  remove(node, false);

  // Add at new place
  Node* current = NULL;
  if ( node->height() < new_height ) {
    if ( node->is_first() ) current = NULL;
    else current = node->previous();
  } else {
    current = first();
  }

  node->set_height(new_height);
  this->add(node, current);
  assert( this->sorted() );
};


// Inefficient, but currently only used for debugging...
size_t NodeContainer::size() const {
  Node* current = first();
  if ( current == NULL ) return 0;

  size_t i = 0;
  for (ConstNodeIterator it = this->iterator(); it.good(); ++it) i++;
  return(i);
};


// Removes all nodes;
// The loop deletes the node from the previous iteration because we still need
// the current node for calling ++it.
void NodeContainer::clear() {
  dout<<"NodeContainer::clear() is called"<<std::endl;
  Node* tmp = NULL;
  for ( NodeIterator it = this->iterator(); it.good(); ++it ) {
    if (tmp != NULL) delete tmp;
    tmp = *it;
  }
  dout<<"all node deleted"<<std::endl;
  if (tmp != NULL) delete tmp;
  set_first(NULL);
  set_last(NULL);
  dout<<"NodeContainer is deleted"<<std::endl;
}

NodeContainer::~NodeContainer() {
  dout<<"first node is "<<this->get(1)<<std::endl;
  dout<<"NodeContainer destructor is called"<<std::endl;
  Node* tmp = NULL;
  for ( NodeIterator it = this->iterator(); it.good(); ++it ) {
    if (tmp != NULL) delete tmp;
    tmp = *it;
  }
  if (tmp != NULL) delete tmp;
  dout<<"all node deleted"<<std::endl;

  set_first(NULL);
  set_last(NULL);
  dout<<"NodeContainer is deleted"<<std::endl;
}

void NodeContainer::add_before(Node* add, Node* next_node){
  //std::cout << "Adding: " << add << " after " << after << std::endl;
  add->set_next(next_node);
  add->set_previous(next_node->previous());

  if ( add->previous() != NULL ) add->previous()->set_next(add);
  next_node->set_previous(add);
  if ( add->is_last() ) this->set_last(add);
}



/*******************************************************
 * Consistency checking
 *******************************************************/

bool NodeContainer::sorted() const {
  Node* current = first();
  if ( !current->is_first() ) {
    std::cout << "NodeContainer: First Node is not first" << std::endl;
    return 0;
  }

  while ( !current->is_last() ) {
    current = current->next();
    if ( current->height() < current->previous()->height() ) {
      std::cout << "NodeContainer: Nodes not sorted" << std::endl;
      return 0;
    }
    if ( current == current->previous() ) {
      std::cout << "NodeContainer: Fatal loop detected" << std::endl;
      return 0;
    }
  }
  
  if ( !current->is_last() ) {
    std::cout << "NodeContainer: Last Node not last" << std::endl;
    return 0;
  }

  return 1;
}


void NodeContainer::print() const {
  //std::cout << "NodeContainer with " << this->size() << " Nodes" << std::endl;
  for ( ConstNodeIterator it = this->iterator(); it.good(); ++it ) {
    std::cout << *it << ": Prev " << (*it)->previous() 
                     << " Next " <<  (*it)->next() << std::endl;
  }
}

#include "node_container.h"

/*******************************************************
 * Con- & Destructor
 *******************************************************/

NodeContainer::NodeContainer() {
  //nodes_.reserve(100);
  set_first(NULL);
  set_last(NULL);
  unsorted_node_ = NULL;
};


NodeContainer::~NodeContainer() {

  if (unsorted_node_ != NULL) delete unsorted_node_;
}



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

Node* NodeContainer::get(size_t nr, const bool &from_sorted) {
  Node* current = first();

  for (size_t i=0; i < nr; ++i) {
    assert(current != NULL);
    current = current->next();
  }
  
  return current;
}


void NodeContainer::add(Node *node, const bool &sort) {
  //std::cout << "First " << first() << std::endl;
  if (first() == NULL) {
    this->set_first(node);
    this->set_last(node);
    return;
  }
  assert(first() != NULL);

  Node* current = first();
  if (node->height() <= current->height()) {
    this->set_first(node);
    node->set_next(current);
    current->set_previous(node);
    return;
  }

  while (node->height() > current->height()) {
    if (current->next() == NULL) break;
    current = current->next();
  }
  
  this->add_after(node, current);
  assert( this->sorted() );
};


void NodeContainer::remove(Node *node) {
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
  assert( this->sorted() );
}


void NodeContainer::move(Node *node, const double &new_height) {
  //Optimize!
  this->remove(node);
  node->set_height(new_height);
  this->add(node);
  assert( this->sorted() );
};


size_t NodeContainer::size() const {
  Node* current = first();
  size_t i = 1;
  while (current->next() != NULL) {
    i++;
    current = current->next();
  }
  return(i);
};


void NodeContainer::add_after(Node* add, Node* after){
  //std::cout << "Adding: " << add << " after " << after << std::endl;
  add->set_previous(after);
  add->set_next(after->next());

  if ( after->next() != NULL ) after->next()->set_previous(add);
  if ( add->is_last() ) this->set_last(add);
  after->set_next(add);
}



/*******************************************************
 * Consistency checking
 *******************************************************/

bool NodeContainer::sorted() const {
  Node* current = first();
  if ( !current->is_first() ) return 0;

  while ( !current->is_last() ) {
    current = current->next();
    if ( current->height() < current->previous()->height() ) return 0;
  }
  
  if ( !current->is_last() ) return 0;
  return 1;
};

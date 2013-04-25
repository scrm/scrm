/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/forest.h"
#include "../src/random/constant_generator.h"

class TestForest : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestForest );
  CPPUNIT_TEST( testInitialization );
  CPPUNIT_TEST( testGettersAndSetters );
  CPPUNIT_TEST( testCreateExampleTree );
  CPPUNIT_TEST( testCheckTreeLength );
  CPPUNIT_TEST( testGetFirstNode );
  CPPUNIT_TEST( testSamplePoint );
  CPPUNIT_TEST( testIsPrunable );
  CPPUNIT_TEST( testPrune );
  //CPPUNIT_TEST( testGetNodeState );
  CPPUNIT_TEST( testPrintTree );
  CPPUNIT_TEST_SUITE_END();

 private:
  Forest *forest;
  ConstantGenerator *rg;

 public:
  void setUp() {
    rg = new ConstantGenerator;
    forest = new Forest(new Model(0), rg);
    forest->createExampleTree();
  }

  void tearDown() {
    delete forest;
  }

  void testInitialization() {
    Forest test_forest = Forest(new Model(4), rg);
    CPPUNIT_ASSERT( test_forest.model().sample_size() == 4 );
    CPPUNIT_ASSERT( test_forest.random_generator() == rg );
    //CPPUNIT_ASSERT( test_forest.local_tree_length() == 0 );
    //CPPUNIT_ASSERT( test_forest.total_tree_length() == 0 );
  }

  void testGettersAndSetters() {
    CPPUNIT_ASSERT( forest->model().sample_size() == 0 );
  }

  void testGetFirstNode() {
    CPPUNIT_ASSERT( forest->nodes()->get(0)->height() == 0 );
  }

  void testCreateExampleTree() {
    CPPUNIT_ASSERT( forest->nodes()->size() == 9 );
    CPPUNIT_ASSERT( forest->local_root() == forest->nodes()->get(8) );
    CPPUNIT_ASSERT( forest->primary_root() == forest->nodes()->get(8) );
    CPPUNIT_ASSERT( forest->local_tree_length() == 24 );
    CPPUNIT_ASSERT( forest->checkTree() == 1 );
  }
  
  void testCheckTreeLength() {
    CPPUNIT_ASSERT( forest->checkTreeLength() );
  }


  void testSamplePoint() {
    TreePoint tp = forest->samplePoint();
    CPPUNIT_ASSERT( tp.base_node() != NULL );
    CPPUNIT_ASSERT( tp.relative_height() > 0 );
  }

  
  void testGetNodeState() { 
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(0), 5) == 0 );
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(8), 5) == 0 );
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(8), 11) == 1 );
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(6), 11) == 2 );
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(7), 11) == 1 );
  }

  void testPrintTree() {
    CPPUNIT_ASSERT_NO_THROW( forest->printTree() );
  }

  void testIsPrunable() {
    forest->writable_model()->set_exact_window_length(5);
    forest->set_current_base(15);
    // Node 6 and 7 are old in this setting, but only node 6 is external and can
    // be pruned.
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(0)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(1)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(2)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(3)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(4)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(5)) );
    CPPUNIT_ASSERT(  forest->isPrunable(forest->getNodes()->get(6)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(7)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(8)) );

    // Orphaned nodes should be pruned
    Node* orphaned = new Node(12, false);
    forest->nodes()->add(orphaned);
    CPPUNIT_ASSERT( forest->isPrunable(orphaned) );

    // Orphaned active not should not occur, but we wont prune them for safety
    // in case they do...
    // Node* orphaned2 = new Node(14, true);
    // forest->nodes()->add(orphaned2);
    // CPPUNIT_ASSERT( !forest->isPrunable(orphaned2) );

    // In-Between Nodes should be pruned, iff they are of same age
    Node *parent = new Node(20, false, 15), 
         *inbetween1 = new Node(19, false, 15), 
         *inbetween2 = new Node(18, false, 13), 
         *child = new Node(17, false, 13);
    forest->nodes()->add(parent);
    forest->nodes()->add(inbetween1);
    forest->nodes()->add(inbetween2);
    forest->nodes()->add(child);
    
    parent->set_first_child(inbetween1);
    inbetween1->set_first_child(inbetween2);
    inbetween2->set_first_child(child);
    child->set_parent(inbetween2);
    inbetween2->set_parent(inbetween1);
    inbetween1->set_parent(parent);
    
    CPPUNIT_ASSERT( !forest->isPrunable(parent) );
    CPPUNIT_ASSERT( !forest->isPrunable(inbetween1) );
    CPPUNIT_ASSERT(  forest->isPrunable(inbetween2) );
    CPPUNIT_ASSERT( !forest->isPrunable(child) );

    // Local in-between nodes
    parent->make_local();
    inbetween1->make_local();
    inbetween2->make_local();
    child->make_local();
    CPPUNIT_ASSERT( !forest->isPrunable(parent) );
    CPPUNIT_ASSERT(  forest->isPrunable(inbetween1) );
    CPPUNIT_ASSERT(  forest->isPrunable(inbetween2) );
    CPPUNIT_ASSERT( !forest->isPrunable(child) );
  }

  void testPrune() {
    // Old node
    forest->writable_model()->set_exact_window_length(5);
    forest->set_current_base(15);
    forest->prune( forest->nodes()->at(6) );
    CPPUNIT_ASSERT( forest->nodes()->size() == 8);
    CPPUNIT_ASSERT( forest->checkTree() );
    // Orphaned node
    forest->prune( forest->nodes()->at(6) );
    CPPUNIT_ASSERT( forest->nodes()->size() == 7);
    CPPUNIT_ASSERT( forest->checkTree() == 1 );

    // In-Between Nodes should be pruned, iff they are of same age
    Node *parent = new Node(20, false, 15), 
         *inbetween1 = new Node(19, false, 15), 
         *inbetween2 = new Node(18, false, 13), 
         *child = new Node(17, false, 13);
    forest->nodes()->add(parent);
    forest->nodes()->add(inbetween1);
    forest->nodes()->add(inbetween2);
    forest->nodes()->add(child);
    
    parent->set_first_child(inbetween1);
    inbetween1->set_first_child(inbetween2);
    inbetween2->set_first_child(child);
    child->set_parent(inbetween2);
    inbetween2->set_parent(inbetween1);
    inbetween1->set_parent(parent);

    forest->prune(inbetween2);
    CPPUNIT_ASSERT( child->parent() == inbetween1 );
    CPPUNIT_ASSERT( inbetween1->parent() == parent );
    CPPUNIT_ASSERT( parent->is_root() );
    CPPUNIT_ASSERT( parent->first_child() == inbetween1 );
    CPPUNIT_ASSERT( inbetween1->first_child() == child );
    CPPUNIT_ASSERT( child->numberOfChildren() == 0 );
  }

};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestForest );

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
  CPPUNIT_TEST( testCheckNodesSorted );
  CPPUNIT_TEST( testGetFirstNode );
  CPPUNIT_TEST( testSamplePoint );
  CPPUNIT_TEST( testCreateRoots );
  CPPUNIT_TEST( testBuildInitialTree );
  CPPUNIT_TEST_SUITE_END();

 private:
  Forest *forest;
  ConstantGenerator *rg;

 public:
  void setUp() {
    rg = new ConstantGenerator;
    forest = new Forest(Model(0), rg);
    forest->createExampleTree();
  }

  void tearDown() {
    delete forest;
  }

  void testInitialization() {
    Forest test_forest = Forest(Model(4), rg);
    CPPUNIT_ASSERT( test_forest.model().sample_size() == 4 );
    CPPUNIT_ASSERT( test_forest.random_generator() == rg );
    CPPUNIT_ASSERT( test_forest.ultimate_root() == NULL );
    CPPUNIT_ASSERT( test_forest.local_tree_length() == 0 );
    CPPUNIT_ASSERT( test_forest.total_tree_length() == 0 );
  }

  void testGettersAndSetters() {
    CPPUNIT_ASSERT( forest->model().sample_size() == 0 );
  }

  void testGetFirstNode() {
    CPPUNIT_ASSERT( forest->getFirstNode()->height() == 0 );
  }

  void testCreateExampleTree() {
    CPPUNIT_ASSERT( forest->countNodes() == 8 );
    CPPUNIT_ASSERT( forest->total_tree_length() == 24 );
    CPPUNIT_ASSERT( forest->checkTree() == 1 );
  }
  
  void testCheckTreeLength() {
    CPPUNIT_ASSERT( forest->checkTreeLength() );
  }

  void testCheckNodesSorted() {
    CPPUNIT_ASSERT(forest->checkNodesSorted() == 1);
    Forest test_forest = Forest(Model(0), rg);
    test_forest.nodes_.push_back(new Node(2));
    test_forest.nodes_.push_back(new Node(1));
    CPPUNIT_ASSERT(test_forest.checkNodesSorted() == 0);
  }

  void testSamplePoint() {
    TreePoint tp = forest->samplePoint();
    CPPUNIT_ASSERT( tp.relative_height() == 4 );
  }

  void testCreateRoots() {
    Forest test_forest = Forest(Model(3), rg);
    test_forest.createRoots();
    CPPUNIT_ASSERT( test_forest.ultimate_root() != NULL );
    CPPUNIT_ASSERT( test_forest.local_root() != NULL );
    Node* u_root = test_forest.ultimate_root();
    Node* l_root = test_forest.local_root();
    CPPUNIT_ASSERT( u_root->height() == FLT_MAX );
    CPPUNIT_ASSERT( l_root->parent() == u_root );
    CPPUNIT_ASSERT( u_root->lower_child() == l_root );
  }
  
  void testBuildInitialTree() {
    //Forest test_forest = Forest(Model(3), rg);
    //test_forest.buildInitialTree();
    //CPPUNIT_ASSERT(test_forest.getFirstNode()->parent() != NULL);
    //CPPUNIT_ASSERT(test_forest.getFirstNode()->parent()->parent() != NULL);
    //CPPUNIT_ASSERT_NO_THROW(forest->checkTree());
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestForest );

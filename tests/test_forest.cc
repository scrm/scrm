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
  CPPUNIT_TEST( testCreateRoots );
  CPPUNIT_TEST( testBuildInitialTree );
  CPPUNIT_TEST( testCountLinesLeft );
  CPPUNIT_TEST( testCountLinesRight );
  CPPUNIT_TEST( testCountBelowLinesLeft );
  CPPUNIT_TEST( testCountBelowLinesRight );
  CPPUNIT_TEST( testPrintTree );
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
    CPPUNIT_ASSERT( forest->nodes()->size() == 8 );
    CPPUNIT_ASSERT( forest->local_root() == forest->nodes()->get(6) );
    CPPUNIT_ASSERT( forest->primary_root() == forest->nodes()->get(6) );
    CPPUNIT_ASSERT( forest->ultimate_root() == forest->nodes()->get(7) );
    //CPPUNIT_ASSERT( forest->total_tree_length() == 24 );
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

  void testCountLinesLeft() {
    CPPUNIT_ASSERT( forest->countLinesLeft(forest->local_root()) == 2 );
    CPPUNIT_ASSERT( forest->countLinesLeft(forest->local_root()->lower_child()) == 1 );
    CPPUNIT_ASSERT( forest->countLinesLeft(forest->nodes()->get(0)) == 0);
  }

  void testCountBelowLinesLeft() {
    CPPUNIT_ASSERT( forest->countBelowLinesLeft(forest->local_root()) == 1 );
    CPPUNIT_ASSERT( forest->countBelowLinesLeft(forest->local_root()->lower_child()) == 0 );
    CPPUNIT_ASSERT( forest->countBelowLinesLeft(forest->nodes()->get(0)) == 0);
  }

  void testCountLinesRight() {
    CPPUNIT_ASSERT( forest->countLinesRight(forest->local_root()) == 2 );
    CPPUNIT_ASSERT( forest->countLinesRight(forest->local_root()->lower_child()) == 1 );
    CPPUNIT_ASSERT( forest->countLinesRight(forest->nodes()->get(0)) == 0);
  }

  void testCountBelowLinesRight() {
    CPPUNIT_ASSERT( forest->countBelowLinesRight(forest->local_root()) == 1 );
    CPPUNIT_ASSERT( forest->countBelowLinesRight(forest->local_root()->lower_child()) == 0 );
    CPPUNIT_ASSERT( forest->countBelowLinesRight(forest->nodes()->get(0)) == 0);
  }

  void testPrintTree() {
    CPPUNIT_ASSERT_NO_THROW( forest->printTree() );
  }

};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestForest );

/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/forest.h"
#include "../src/node.h"
#include "../src/random/constant_generator.h"

class TestForest : public CppUnit::TestCase {

	CPPUNIT_TEST_SUITE( TestForest );
	CPPUNIT_TEST( testInitialization );
 	CPPUNIT_TEST( testGettersAndSetters );
 	CPPUNIT_TEST( testGetFirstNode );
 	CPPUNIT_TEST( testBuildInitialTree );
	CPPUNIT_TEST_SUITE_END();

    private:
     Node *node1, *node2, *node3, *node4;
     Forest *forest;

	public:
      void setUp() {
         node1 = new Node(1);
         node2 = new Node(2);
         node3 = new Node(3);
         node4 = new Node(4);
         
         forest = new Forest(Model(0), new ConstantGenerator());
         
         forest->addNode(node1);
         forest->addNode(node2);
         forest->addNode(node3);
         forest->addNode(node4);
      }

      void tearDown() {
         delete node1;
         delete node2; 
         delete node3; 
         delete node4;
         delete forest;
      }

      void testInitialization() {
          CPPUNIT_ASSERT( forest->countNodes() == 4 );
      }

      void testGettersAndSetters() {
          CPPUNIT_ASSERT( forest->model().sample_size() == 0 );
      }

      void testGetFirstNode() {
          CPPUNIT_ASSERT( forest->getFirstNode()->height() == 1 );
      }

      void testBuildInitialTree() {
          Forest test_forest = Forest(Model(3), new ConstantGenerator());
          test_forest.buildInitialTree();
          CPPUNIT_ASSERT(test_forest.getFirstNode()->parent() != NULL);
          CPPUNIT_ASSERT(test_forest.getFirstNode()->parent()->parent() != NULL);
      }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestForest );

/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/forest.h"
#include "../src/node.h"

class TestForest : public CppUnit::TestCase {

	CPPUNIT_TEST_SUITE( TestForest );

	CPPUNIT_TEST( testAddNode );
 	CPPUNIT_TEST( testGettersAndSetters );
 	CPPUNIT_TEST( testCreateSampleNodes );
	CPPUNIT_TEST_SUITE_END();

	public:
        friend class Forest;
      
		void testAddNode() {
            Forest forest = Forest();
            Node* node1 = new Node();
            Node* node2 = new Node();
            forest.addNode(node1);
			CPPUNIT_ASSERT( forest.countNodes() == 1 );
            forest.addNode(node2);
			CPPUNIT_ASSERT( forest.countNodes() == 2 );
		}

        void testGettersAndSetters() {
            Forest forest = Forest();
            forest.set_model(Model(7));
			CPPUNIT_ASSERT( forest.model().sample_size() == 7 );
        }

        void testCreateSampleNodes() {
            Forest forest = Forest(Model(10));
            forest.createSampleNodes();
			CPPUNIT_ASSERT( forest.countNodes() == 10 );
        }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestForest );

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

	CPPUNIT_TEST_SUITE_END();

	public:
		void testAddNode() {
            Forest forest = Forest();
            Node node1, node2;
            forest.addNode(node1);
			CPPUNIT_ASSERT( forest.countNodes() == 1 );
            forest.addNode(node2);
			CPPUNIT_ASSERT( forest.countNodes() == 2 );
		}
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestForest );

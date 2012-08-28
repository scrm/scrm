/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/node.h"

class TestNode : public CppUnit::TestCase {

	CPPUNIT_TEST_SUITE( TestNode );

	CPPUNIT_TEST( testGettersAndSetters );

	CPPUNIT_TEST_SUITE_END();

	public:
		void testGettersAndSetters() {
            Node node1, node2;

            //height
            node1.set_height(1);
			CPPUNIT_ASSERT( node1.get_height() == 1 );

			//parent
            node2.set_parent(&node1);
			CPPUNIT_ASSERT( node2.get_parent()->get_height() == 1 );

            //Children
            node2.set_higher_child(&node1);
			CPPUNIT_ASSERT( node2.get_higher_child()->get_height() == 1 );
            node2.set_lower_child(&node1);
			CPPUNIT_ASSERT( node2.get_lower_child()->get_height() == 1 );
		}
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestNode );

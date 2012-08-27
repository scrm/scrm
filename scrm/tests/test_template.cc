/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>


class MyTest : public CppUnit::TestCase {

	CPPUNIT_TEST_SUITE( MyTest );

	CPPUNIT_TEST( testOneThing );
	CPPUNIT_TEST( testAnother );
	CPPUNIT_TEST( testException );

	CPPUNIT_TEST_SUITE_END();

	public:
		void testOneThing() {
			CPPUNIT_ASSERT( 1 == 1 );
			CPPUNIT_ASSERT_MESSAGE("error message", 1 == 1 );
		}
		
		void testAnother() {
			CPPUNIT_ASSERT_EQUAL(1, 1);
		}

		void testException() {
			CPPUNIT_ASSERT( 1 == 1 );
			CPPUNIT_ASSERT_THROW(
				throw std::exception(), std::exception );
		}
};

CPPUNIT_TEST_SUITE_REGISTRATION( MyTest );

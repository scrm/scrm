/*
 * A sample test case which can be used as a template.
 */
#include <iostream>
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/model.h"

class TestModel : public CppUnit::TestCase {

	CPPUNIT_TEST_SUITE( TestModel );

 	CPPUNIT_TEST( testGettersAndSetters );
	CPPUNIT_TEST_SUITE_END();

	public:
        void testGettersAndSetters() {
            Model model = Model();
            
            model.setSampleSize(5);
			CPPUNIT_ASSERT( model.getSampleSize() == 5 );

            model.setPopulationSize(5000);
			CPPUNIT_ASSERT( model.getPopulationSize() == 5000 );

            model.setMutationRate(10.7);
			CPPUNIT_ASSERT( model.getMutationRate() == 10.7 );

            model.setRecombinationRate(20.1);
			CPPUNIT_ASSERT( model.getRecombinationRate() == 20.1 );
        }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestModel );

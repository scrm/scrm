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
            
            model.set_sample_size(5);
			CPPUNIT_ASSERT( model.sample_size() == 5 );

            model.set_population_size(5000);
			CPPUNIT_ASSERT( model.population_size() == 5000 );

            model.set_mutation_rate(10.7);
			CPPUNIT_ASSERT( model.mutation_rate() == 10.7 );

            model.set_recombination_rate(20.1);
			CPPUNIT_ASSERT( model.recombination_rate() == 20.1 );
        }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestModel );

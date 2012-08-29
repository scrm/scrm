/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/random.h"

class TestRandomGenerator : public CppUnit::TestCase {

	CPPUNIT_TEST_SUITE( TestRandomGenerator );

	CPPUNIT_TEST( testSeeding );
	CPPUNIT_TEST( testSampleInt );
	CPPUNIT_TEST( testSampleTwoElements );
	CPPUNIT_TEST_SUITE_END();

	public:
	  void testSeeding() {
            RandomGenerator rg = RandomGenerator();
            rg.set_seed(123);
            double sample1 = rg.sample();
            rg.set_seed(123);
            double sample2 = rg.sample();
			CPPUNIT_ASSERT( sample1 == sample2 );
		}

      //Only very week tests below, mostly to assure that the functions actually
      //return something...
      void testSample() {
        RandomGenerator rg = RandomGenerator();
        double sample = rg.sample();
		CPPUNIT_ASSERT( 0 <= sample );
		CPPUNIT_ASSERT( sample < 1 );
      }

      void testSampleInt() {
        RandomGenerator rg = RandomGenerator();
        int sample = rg.sampleInt(10);
		CPPUNIT_ASSERT( 0 <= sample );
		CPPUNIT_ASSERT( sample < 10 );
      }

      void testSampleTwoElements() {
        RandomGenerator rg = RandomGenerator();
        int sample1, sample2;
        rg.sampleTwoElements(2, &sample1, &sample2);
		CPPUNIT_ASSERT( sample1 != sample2 );
      }

};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestRandomGenerator );

/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/random/mersenne_twister.h"
#include "../src/random/random_generator.h"

class TestRandomGenerator : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestRandomGenerator );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testSampleExpoLimit );

  CPPUNIT_TEST_SUITE_END();

 private:
  RandomGenerator *rg;

 public:
  void setUp() {
    rg = new MersenneTwister(5);
  }

  void tearDown() {
    delete rg;
  }

  void testConstructor() {
    CPPUNIT_ASSERT_EQUAL( (int)5, rg->seed() );
    CPPUNIT_ASSERT( rg->unit_exponential_ > 0 );
  }

  void testSampleExpoLimit() {
    double expo;
    for (size_t i = 0; i < 1000; ++i) {
      expo = rg->sampleExpoLimit(1, 2);
      CPPUNIT_ASSERT( expo == -1 || expo > 0 );
      CPPUNIT_ASSERT( expo <= 2 ); 
    }
  }

};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestRandomGenerator );

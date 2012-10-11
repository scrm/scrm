/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/forest.h"
#include "../src/random/constant_generator.h"


class TestEvent : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestEvent );

  CPPUNIT_TEST( testEventIteratorCreation );

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

  void testEventIteratorCreation() {
    EventIterator events = EventIterator(forest, 0);
    events = EventIterator(forest, 2);
    
    //CPPUNIT_ASSERT( 1 == 1 );
    //CPPUNIT_ASSERT_MESSAGE("error message", 1 == 1 );
  }

};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestEvent );

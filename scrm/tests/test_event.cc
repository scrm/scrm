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
  CPPUNIT_TEST( testEventIteratorNext );

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
    EventIterator it = EventIterator(forest, 0);
    CPPUNIT_ASSERT( it.start_height_ == 0 );
    CPPUNIT_ASSERT( it.contemporaries_.size() == 0 );
    CPPUNIT_ASSERT( *it.twig_iterator_ == forest->getFirstNode() );
    it = EventIterator(forest, 2);
    CPPUNIT_ASSERT( it.start_height_ == 2 );
    CPPUNIT_ASSERT( it.contemporaries_.size() == 3 );
    it = EventIterator(forest, 4);
    CPPUNIT_ASSERT( it.start_height_ == 4 );
    CPPUNIT_ASSERT( it.contemporaries_.size() == 2 );
  }

  void testEventIteratorNext() {
    EventIterator events = EventIterator(forest, 0.5);
    Event event = events.next();
    CPPUNIT_ASSERT( event.start_height() == 0.5 );
    CPPUNIT_ASSERT( event.end_height() == 1 );
    CPPUNIT_ASSERT( event.contemporaries().size() == 4);
    event = events.next();
    CPPUNIT_ASSERT( event.start_height() == 1 );
    CPPUNIT_ASSERT( event.end_height() == 3 );
    CPPUNIT_ASSERT( event.contemporaries().size() == 3);
    event = events.next();
    CPPUNIT_ASSERT( event.start_height() == 3 );
    CPPUNIT_ASSERT( event.end_height() == 10 );
    CPPUNIT_ASSERT( event.contemporaries().size() == 2);
    event = events.next();
    CPPUNIT_ASSERT( event.start_height() == 10 );
    CPPUNIT_ASSERT( event.end_height() == FLT_MAX );
    CPPUNIT_ASSERT( event.contemporaries().size() == 1);
    CPPUNIT_ASSERT_THROW( event = events.next(), std::out_of_range );

    events = EventIterator(forest, 0);
    event = events.next();
    CPPUNIT_ASSERT( event.start_height() == 0 );
    CPPUNIT_ASSERT( event.end_height() == 1 );
    CPPUNIT_ASSERT( event.contemporaries().size() == 4 );
  }

};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestEvent );

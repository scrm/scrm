/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/forest.h"
#include "../src/random/constant_generator.h"


class TestTimeInterval : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestTimeInterval );

  CPPUNIT_TEST( testIteratorCreation );
  CPPUNIT_TEST( testIteratorNext );
  CPPUNIT_TEST( testSplitIntervall );
  CPPUNIT_TEST( testSampleContemporary );
  CPPUNIT_TEST( testRecalculateTI );

  CPPUNIT_TEST_SUITE_END();

 private:
  Forest *forest;
  ConstantGenerator *rg;

 public:
  void setUp() {
    rg = new ConstantGenerator();
    forest = new Forest(Model(0), new MersenneTwister(5));
    forest->createExampleTree();
  }

  void tearDown() {
    delete forest;
  }

  void testIteratorCreation() {
    TimeIntervalIterator it = TimeIntervalIterator(forest, forest->getNodes()->get(0));
    CPPUNIT_ASSERT( (*it).start_height() == 0 );
    CPPUNIT_ASSERT( (*it).end_height() == 1 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 4 );
    
    it = TimeIntervalIterator(forest, forest->getNodes()->get(4));
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 3 );
    
    it = TimeIntervalIterator(forest, forest->getNodes()->get(5));
    CPPUNIT_ASSERT( (*it).start_height() == 3 );
    CPPUNIT_ASSERT( (*it).end_height() == 4 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 2 );
 
    it = TimeIntervalIterator(forest, forest->getNodes()->get(6));
    CPPUNIT_ASSERT( (*it).start_height() == 4 );
    CPPUNIT_ASSERT( (*it).end_height() == 6 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 3 );
 
    it = TimeIntervalIterator(forest, forest->getNodes()->get(7));
    CPPUNIT_ASSERT( (*it).start_height() == 6 );
    CPPUNIT_ASSERT( (*it).end_height() == 10 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 2 );
    
    it = TimeIntervalIterator(forest, forest->getNodes()->get(8));
    CPPUNIT_ASSERT( (*it).start_height() == 10 );
    CPPUNIT_ASSERT( (*it).end_height() == FLT_MAX );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 0 );
  }

  void testIteratorNext() {
    TimeIntervalIterator it = TimeIntervalIterator(forest, forest->getNodes()->get(0));
    CPPUNIT_ASSERT( (*it).start_height() == 0 );
    CPPUNIT_ASSERT( (*it).end_height() == 1 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 4);
    CPPUNIT_ASSERT( it.good() );
    
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 3);
    CPPUNIT_ASSERT( it.good() );
  
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 3 );
    CPPUNIT_ASSERT( (*it).end_height() == 4 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 2 );
    CPPUNIT_ASSERT( it.good() );
 
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 4 );
    CPPUNIT_ASSERT( (*it).end_height() == 6 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 3 );
    CPPUNIT_ASSERT( it.good() );
 
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 6 );
    CPPUNIT_ASSERT( (*it).end_height() == 10 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 2 );
    CPPUNIT_ASSERT( it.good() );
    
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 10 );
    CPPUNIT_ASSERT( (*it).end_height() == FLT_MAX );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 0 );
    CPPUNIT_ASSERT( it.good() );
  
    ++it;
    CPPUNIT_ASSERT( !it.good() );

    int i = 0;
    for (TimeIntervalIterator events = TimeIntervalIterator(forest, forest->getNodes()->get(0));
         events.good(); ++events) { 
      ++i;
      //std::cout << (*events).start_height() << " - " << (*events).end_height() << std::endl; 
    }
    CPPUNIT_ASSERT( i == 6 );
  }

  void testSplitIntervall() {
    TimeIntervalIterator it = TimeIntervalIterator(forest, forest->getNodes()->get(0));
    Node* split_node = new Node(0.5);
    it.splitCurrentInterval(split_node);
    ++it;

    // Check split interval
    CPPUNIT_ASSERT( (*it).start_height() == 0.5 );
    CPPUNIT_ASSERT( (*it).end_height() == 1 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 4);
    CPPUNIT_ASSERT( it.good() );

    // Check that next interval is unchanged
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 3);
    CPPUNIT_ASSERT( it.good() );
  }

  void testSampleContemporary() {
    TimeIntervalIterator it = TimeIntervalIterator(forest, forest->getNodes()->get(7));
    CPPUNIT_ASSERT( (*it).getRandomContemporary() != NULL );
    CPPUNIT_ASSERT( (*it).getRandomContemporary() != NULL );
    CPPUNIT_ASSERT( (*it).getRandomContemporary() != NULL );
    CPPUNIT_ASSERT( (*it).getRandomContemporary() != NULL );
    CPPUNIT_ASSERT( (*it).getRandomContemporary() != NULL );
    CPPUNIT_ASSERT( (*it).getRandomContemporary() != NULL );
    CPPUNIT_ASSERT( (*it).getRandomContemporary() != NULL );
    CPPUNIT_ASSERT( (*it).getRandomContemporary() != NULL );
    CPPUNIT_ASSERT( (*it).getRandomContemporary() != NULL );
    CPPUNIT_ASSERT( (*it).getRandomContemporary() != NULL );

    ++it;
    CPPUNIT_ASSERT_THROW( (*it).getRandomContemporary(), std::out_of_range );
  }

  void testRecalculateTI() {
    TimeIntervalIterator it = TimeIntervalIterator(forest, forest->getNodes()->get(0));
    ++it;

    // Check that is does not change anything if the underlying tree is
    // unchanged
    it.recalculateInterval();
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 3);
    CPPUNIT_ASSERT( it.good() );

    // Now change the current interval
    forest->nodes()->add(new Node(2));
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    it.recalculateInterval();
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 2 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 3);
    CPPUNIT_ASSERT( it.good() );

    // Also test the next interval
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 2 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( (*it).numberOfContemporaries() == 3);
    CPPUNIT_ASSERT( it.good() );
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestTimeInterval );

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
  CPPUNIT_TEST( testSplitIntervall );

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
    EventIterator it = EventIterator(forest, forest->getNodes()->get(0));
    CPPUNIT_ASSERT( (*it).start_height() == 0 );
    CPPUNIT_ASSERT( (*it).end_height() == 1 );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 4 );
    
    it = EventIterator(forest, forest->getNodes()->get(4));
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 3 );
    
    it = EventIterator(forest, forest->getNodes()->get(5));
    CPPUNIT_ASSERT( (*it).start_height() == 3 );
    CPPUNIT_ASSERT( (*it).end_height() == 4 );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 2 );
 
    it = EventIterator(forest, forest->getNodes()->get(6));
    CPPUNIT_ASSERT( (*it).start_height() == 4 );
    CPPUNIT_ASSERT( (*it).end_height() == 6 );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 3 );
 
    it = EventIterator(forest, forest->getNodes()->get(7));
    CPPUNIT_ASSERT( (*it).start_height() == 6 );
    CPPUNIT_ASSERT( (*it).end_height() == 10 );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 2 );
    
    it = EventIterator(forest, forest->getNodes()->get(8));
    CPPUNIT_ASSERT( (*it).start_height() == 10 );
    CPPUNIT_ASSERT( (*it).end_height() == FLT_MAX );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 0 );
  }

  void testEventIteratorNext() {
    EventIterator it = EventIterator(forest, forest->getNodes()->get(0));
    CPPUNIT_ASSERT( (*it).start_height() == 0 );
    CPPUNIT_ASSERT( (*it).end_height() == 1 );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 4);
    CPPUNIT_ASSERT( it.good() );
    
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 3);
    CPPUNIT_ASSERT( it.good() );
  
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 3 );
    CPPUNIT_ASSERT( (*it).end_height() == 4 );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 2 );
    CPPUNIT_ASSERT( it.good() );
 
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 4 );
    CPPUNIT_ASSERT( (*it).end_height() == 6 );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 3 );
    CPPUNIT_ASSERT( it.good() );
 
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 6 );
    CPPUNIT_ASSERT( (*it).end_height() == 10 );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 2 );
    CPPUNIT_ASSERT( it.good() );
    
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 10 );
    CPPUNIT_ASSERT( (*it).end_height() == FLT_MAX );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 0 );
    CPPUNIT_ASSERT( it.good() );
  
    ++it;
    CPPUNIT_ASSERT( !it.good() );

    int i = 0;
    for (EventIterator events = EventIterator(forest, forest->getNodes()->get(0));
         events.good(); ++events) { 
      ++i;
      //std::cout << (*events).start_height() << " - " << (*events).end_height() << std::endl; 
    }
    CPPUNIT_ASSERT( i == 6 );
  }

  void testSplitIntervall() {
    EventIterator it = EventIterator(forest, forest->getNodes()->get(0));
    Node* split_node = new Node(0.5);
    it.splitCurrentInterval(split_node);
    ++it;

    // Check split interval
    CPPUNIT_ASSERT( (*it).start_height() == 0.5 );
    CPPUNIT_ASSERT( (*it).end_height() == 1 );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 4);
    CPPUNIT_ASSERT( it.good() );

    // Check that next interval is unchanged
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( (*it).contemporaries().size() == 3);
    CPPUNIT_ASSERT( it.good() );
  }

};

//Uncomment this to activate the test
// CPPUNIT_TEST_SUITE_REGISTRATION( TestEvent );

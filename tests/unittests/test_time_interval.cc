/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../../src/forest.h"
#include "../../src/random/mersenne_twister.h"


class TestTimeInterval : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestTimeInterval );

  CPPUNIT_TEST( testIteratorCreation );
  CPPUNIT_TEST( testIteratorNext );
  CPPUNIT_TEST( testIteratorCreationWithTimeFrames );
  CPPUNIT_TEST( testIteratorNextWithTimeFrames );
  CPPUNIT_TEST( testSplitIntervall );
  CPPUNIT_TEST( testRecalculateTI );
  CPPUNIT_TEST( testFindContemporaries );

  CPPUNIT_TEST_SUITE_END();

 private:
  Forest *forest, *forest2;
  Model *model, *model2;
  MersenneTwister *rg;

 public:
  void setUp() {
    rg = new MersenneTwister(5);
    model = new Model(5);
    model2 = new Model(5);
    model2->set_population_number(3);

    forest = new Forest(model, rg);
    forest->createExampleTree();

    forest2 = new Forest(model2, rg);
    forest2->createExampleTree();
  }

  void tearDown() {
    delete forest, forest2;
    delete model, model2; 
    delete rg;
  }

  void testIteratorCreation() {
    TimeIntervalIterator it(forest, forest->nodes()->at(0));
    CPPUNIT_ASSERT( (*it).start_height() == 0 );
    CPPUNIT_ASSERT( (*it).end_height() == 1 );
    CPPUNIT_ASSERT_EQUAL( (size_t)4, it.contemporaries()->size(0) );
    
    TimeIntervalIterator it2(forest, forest->nodes()->at(4));
    CPPUNIT_ASSERT( (*it2).start_height() == 1 );
    CPPUNIT_ASSERT( (*it2).end_height() == 3 );
    CPPUNIT_ASSERT_EQUAL( (size_t)3, it2.contemporaries()->size(0) );
    
    TimeIntervalIterator it3(forest, forest->nodes()->at(5));
    CPPUNIT_ASSERT( (*it3).start_height() == 3 );
    CPPUNIT_ASSERT( (*it3).end_height() == 4 );
    CPPUNIT_ASSERT( it3.contemporaries()->size(0) == 2 );
 
    TimeIntervalIterator it4(forest, forest->nodes()->at(6));
    CPPUNIT_ASSERT( (*it4).start_height() == 4 );
    CPPUNIT_ASSERT( (*it4).end_height() == 6 );
    CPPUNIT_ASSERT( it4.contemporaries()->size(0) == 3 );
 
    TimeIntervalIterator it5(forest, forest->nodes()->at(7));
    CPPUNIT_ASSERT( (*it5).start_height() == 6 );
    CPPUNIT_ASSERT( (*it5).end_height() == 10 );
    CPPUNIT_ASSERT( it5.contemporaries()->size(0) == 2 );
    
    TimeIntervalIterator it6(forest, forest->nodes()->at(8));
    CPPUNIT_ASSERT( (*it6).start_height() == 10 );
    CPPUNIT_ASSERT( (*it6).end_height() == DBL_MAX );
    CPPUNIT_ASSERT( it6.contemporaries()->size(0) == 0 );
  }

  void testIteratorNext() {
    TimeIntervalIterator it(forest, forest->nodes()->at(0));
    CPPUNIT_ASSERT( (*it).start_height() == 0 );
    CPPUNIT_ASSERT( (*it).end_height() == 1 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 4);
    CPPUNIT_ASSERT( it.good() );
    
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 3);
    CPPUNIT_ASSERT( it.good() );
  
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 3 );
    CPPUNIT_ASSERT( (*it).end_height() == 4 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 2 );
    CPPUNIT_ASSERT( it.good() );
 
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 4 );
    CPPUNIT_ASSERT( (*it).end_height() == 6 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 3 );
    CPPUNIT_ASSERT( it.good() );
 
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 6 );
    CPPUNIT_ASSERT( (*it).end_height() == 10 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 2 );
    CPPUNIT_ASSERT( it.good() );
    
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 10 );
    CPPUNIT_ASSERT( (*it).end_height() == DBL_MAX );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 0 );
    CPPUNIT_ASSERT( it.good() );
  
    ++it;
    CPPUNIT_ASSERT( !it.good() );

    int i = 0;
    for (TimeIntervalIterator events(forest, forest->nodes()->at(0));
         events.good(); ++events) { 
      ++i;
    }
    CPPUNIT_ASSERT( i == 6 );
  }

  void testIteratorCreationWithTimeFrames() {
    forest->writable_model()->addGrowthRates(0.0, std::vector<double>(1, 1.5));
    forest->writable_model()->addGrowthRates(0.5, std::vector<double>(1, 1.5));
    forest->writable_model()->addGrowthRates(1.5, std::vector<double>(1, 1.5));

    TimeIntervalIterator it(forest, forest->nodes()->at(0));
    CPPUNIT_ASSERT( (*it).start_height() == 0 );
    CPPUNIT_ASSERT( (*it).end_height() == 0.5 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 4 );

    TimeIntervalIterator it2(forest, forest->nodes()->at(4));
    CPPUNIT_ASSERT( (*it2).start_height() == 1 );
    CPPUNIT_ASSERT( (*it2).end_height() == 1.5 );
    CPPUNIT_ASSERT_EQUAL( (size_t)3, it2.contemporaries()->size(0) );
  }

  void testIteratorNextWithTimeFrames() {
    forest->writable_model()->addGrowthRates(0.0, std::vector<double>(1, 1.5));
    forest->writable_model()->addGrowthRates(0.5, std::vector<double>(1, 1.5));
    forest->writable_model()->addGrowthRates(1, std::vector<double>(1, 1.5));
    forest->writable_model()->addGrowthRates(1.5, std::vector<double>(1, 1.5));

    TimeIntervalIterator it(forest, forest->nodes()->at(0));
    CPPUNIT_ASSERT( (*it).start_height() == 0 );
    CPPUNIT_ASSERT( (*it).end_height() == 0.5 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 4 );

    CPPUNIT_ASSERT_NO_THROW( it.next() );
    CPPUNIT_ASSERT( (*it).start_height() == 0.5 );
    CPPUNIT_ASSERT( (*it).end_height() == 1 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 4 );

    CPPUNIT_ASSERT_NO_THROW( it.next() );
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 1.5 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 3 );
    
    CPPUNIT_ASSERT_NO_THROW( it.next() );
    CPPUNIT_ASSERT( (*it).start_height() == 1.5 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 3 );
  }

  void testSplitIntervall() {
    TimeIntervalIterator it(forest, forest->nodes()->at(0));
    Node* split_node = new Node(0.5);
    it.splitCurrentInterval(split_node);
    ++it;

    // Check split interval
    CPPUNIT_ASSERT( (*it).start_height() == 0.5 );
    CPPUNIT_ASSERT( (*it).end_height() == 1 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 4);
    CPPUNIT_ASSERT( it.good() );

    // Check that next interval is unchanged
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 3);
    CPPUNIT_ASSERT( it.good() );
    delete split_node;
  }

  void testRecalculateTI() {
    TimeIntervalIterator it(forest, forest->nodes()->at(0));
    ++it;

    // Check that is does not change anything if the underlying tree is
    // unchanged
    it.recalculateInterval();
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 3);
    CPPUNIT_ASSERT( it.good() );

    // Now change the current interval
    forest->nodes()->add(new Node(2));
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    it.recalculateInterval();
    CPPUNIT_ASSERT( (*it).start_height() == 1 );
    CPPUNIT_ASSERT( (*it).end_height() == 2 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 3);
    CPPUNIT_ASSERT( it.good() );

    // Also test the next interval
    ++it;
    CPPUNIT_ASSERT( (*it).start_height() == 2 );
    CPPUNIT_ASSERT( (*it).end_height() == 3 );
    CPPUNIT_ASSERT( it.contemporaries()->size(0) == 3);
    CPPUNIT_ASSERT( it.good() );

    // Test if we can change the final interval
    while ( (*it).end_height() < 1000 ) ++it;
    forest->nodes()->add(new Node(1000));
    it.recalculateInterval();
    CPPUNIT_ASSERT( (*it).end_height() == 1000 ); 
    CPPUNIT_ASSERT_NO_THROW( ++it );
    CPPUNIT_ASSERT( (*it).start_height() == 1000 ); 
    CPPUNIT_ASSERT( (*it).end_height() == DBL_MAX ); 
  }

  void testFindContemporaries() {
    TimeIntervalIterator tii(forest, forest->nodes()->at(0));
    tii.contemporaries()->clear();

    tii.searchContemporariesBottomUp(forest->nodes()->at(4));
    CPPUNIT_ASSERT_EQUAL((size_t)2, tii.contemporaries()->size(0)); 
    tii.contemporaries()->buffer(forest->nodes()->at(4)->height());

    tii.searchContemporariesBottomUp(forest->nodes()->at(5));
    CPPUNIT_ASSERT_EQUAL((size_t)1, tii.contemporaries()->size(0)); 
    tii.searchContemporariesBottomUp(forest->nodes()->at(5), true);
    CPPUNIT_ASSERT_EQUAL((size_t)1, tii.contemporaries()->size(0)); 
    tii.contemporaries()->buffer(forest->nodes()->at(4)->height());

    tii.searchContemporariesBottomUp(forest->nodes()->at(6));
    CPPUNIT_ASSERT_EQUAL((size_t)2, tii.contemporaries()->size(0)); 
    tii.searchContemporariesBottomUp(forest->nodes()->at(6), true);
    CPPUNIT_ASSERT_EQUAL((size_t)2, tii.contemporaries()->size(0)); 

    tii.searchContemporariesBottomUp(forest->nodes()->at(7));
    CPPUNIT_ASSERT_EQUAL((size_t)2, tii.contemporaries()->size(0)); 
    tii.searchContemporariesBottomUp(forest->nodes()->at(7), true);
    CPPUNIT_ASSERT_EQUAL((size_t)2, tii.contemporaries()->size(0)); 

    tii.searchContemporariesBottomUp(forest->nodes()->at(8));
    CPPUNIT_ASSERT_EQUAL((size_t)0, tii.contemporaries()->size(0)); 
    tii.searchContemporariesBottomUp(forest->nodes()->at(8), true);
    CPPUNIT_ASSERT_EQUAL((size_t)0, tii.contemporaries()->size(0)); 
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestTimeInterval );

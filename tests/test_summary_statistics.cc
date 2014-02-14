/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <valarray>

#include "../src/forest.h"
#include "../src/tree_point.h"
#include "../src/random/mersenne_twister.h"
#include "../src/summary_statistics/tmrca.h"
#include "../src/summary_statistics/seg_sites.h"
#include "../src/summary_statistics/summary_statistic.h"

class TestSummaryStatistics : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestSummaryStatistics );

  CPPUNIT_TEST( testTMRCA );
  CPPUNIT_TEST( testSegSitesGetHaplotypes );
  CPPUNIT_TEST( testSegSitesCalculate );

  CPPUNIT_TEST_SUITE_END();

 private:
  Forest *forest;
  MersenneTwister *rg;

 public:
  void setUp() {
    rg = new MersenneTwister(1234);
    forest = new Forest(new Model(5), rg);
    forest->createExampleTree();
  }

  void tearDown() {
    delete forest->writable_model();
    delete forest;
    delete rg;
  }

  void testTMRCA() {
    TMRCA* tmrca = new TMRCA();
    tmrca->calculate(*forest);
    delete tmrca;

    SummaryStatistic* tmrca2 = new TMRCA();
    tmrca2->calculate(*forest);
    delete tmrca2;
  }

  void testSegSitesGetHaplotypes() {
    SegSites seg_sites = SegSites();

    TreePoint mutation = TreePoint(forest->nodes()->at(4), 1, true);
    valarray<bool> haplotype = seg_sites.getHaplotypes(mutation, *forest);
    CPPUNIT_ASSERT_EQUAL( (size_t)4, haplotype.size() );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[3] );

    mutation = TreePoint(forest->nodes()->at(5), 0.5, true);
    haplotype = seg_sites.getHaplotypes(mutation, *forest);
    CPPUNIT_ASSERT_EQUAL( (size_t)4, haplotype.size() );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[3] );

    mutation = TreePoint(forest->nodes()->at(8), 0.5, true);
    haplotype = seg_sites.getHaplotypes(mutation, *forest);
    CPPUNIT_ASSERT_EQUAL( (size_t)4, haplotype.size() );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[3] );
  }

  void testSegSitesCalculate() {
    forest->createScaledExampleTree();
    forest->writable_model()->set_mutation_rate(0.0001);
    SegSites seg_sites = SegSites();
    seg_sites.calculate(*forest);
    CPPUNIT_ASSERT( seg_sites.countMutations() > 0 ); 
    int freqs[4] = { 0 };
    int types[2] = { 0 };
    int sum = 0;
    for (auto it = seg_sites.haplotypes_.begin(); it != seg_sites.haplotypes_.end(); ++it) {
      for (size_t i = 0; i < 4; ++i) {
        if ((*it)[i]) ++freqs[i]; 
        sum += (*it)[i];
      }
      CPPUNIT_ASSERT( sum == 1 || sum == 2 );
      ++types[sum-1];
      sum = 0;
    }

    for (size_t i = 0; i < 4; ++i) {
      CPPUNIT_ASSERT( 0.4 < (double)freqs[i] / seg_sites.countMutations() );
      CPPUNIT_ASSERT( (double)freqs[i] / seg_sites.countMutations() < 0.43 );
    }

    CPPUNIT_ASSERT( (double)types[0] / types[1] > 0.45 );
    CPPUNIT_ASSERT( (double)types[0] / types[1] < 0.55 );
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestSummaryStatistics );

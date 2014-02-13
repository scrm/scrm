/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/forest.h"
#include "../src/random/mersenne_twister.h"
#include "../src/summary_statistics/tmrca.h"
#include "../src/summary_statistics/summary_statistic.h"

class TestSummaryStatistics : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestSummaryStatistics );

  CPPUNIT_TEST( testTMRCA );

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

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestSummaryStatistics );

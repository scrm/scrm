#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../../src/random/fastfunc.h"
#include "../../src/random/mersenne_twister.h"

class TestFastFunc : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestFastFunc );

  CPPUNIT_TEST( testlog );
  CPPUNIT_TEST( testexp );

  CPPUNIT_TEST_SUITE_END();

 public:
  void testlog() {
    MersenneTwister rg(5);
    class FastFunc ff;
    double d = 1e-6;
    double maxdiff = 0.0;

    // Check difference to log
    double logd, truelogd;
    while (d < 1) {
      logd = ff.fastlog( d );
      CPPUNIT_ASSERT( logd != 0.0 );
      truelogd = log( d );
      maxdiff = std::max( maxdiff, std::abs( logd-truelogd ) );
      d += 1e-6;
    }
    CPPUNIT_ASSERT( maxdiff < 1.5e-7 );

    // Check that the random generator does not produce exactly 0
    for (int i=0; i<1e6; ++i) CPPUNIT_ASSERT( -ff.fastlog(rg.sample()) != 0.0 );
  }

  void testexp() {
    class MersenneTwister rg;
    class FastFunc ff;
    int i;

    double x, true_exp, lower_bound, upper_bound;

    for (i=0; i<10000; ++i) {
      x = (rg.sample()-0.5) * 1400.0;
      true_exp = exp(x);
      lower_bound = ff.fastexp_lo(x);
      upper_bound = ff.fastexp_up(x);
      //std::cout << x << " " << true_exp << " " << lower_bound << " " << upper_bound << std::endl;
      CPPUNIT_ASSERT(lower_bound <= true_exp);
      CPPUNIT_ASSERT(lower_bound > true_exp * (1.0 - 0.05792));
      CPPUNIT_ASSERT(upper_bound >= true_exp);
      CPPUNIT_ASSERT(upper_bound < true_exp * (1.0 + 0.06148));
    }
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestFastFunc );

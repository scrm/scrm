#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../src/random/fastfunc.h"
#include "../src/random/mersenne_twister.h"

class TestFastFunc : public CppUnit::TestCase {

	CPPUNIT_TEST_SUITE( TestFastFunc );

	CPPUNIT_TEST( testlog );
	CPPUNIT_TEST( testexp );

	CPPUNIT_TEST_SUITE_END();

	public:
    void testlog() {
      class MersenneTwister rg;
      class FastFunc ff;
      double d = 1.0;
      double maxdiff = 0.0;

      double logd, truelogd;
      for (int i=0; i<10000000; i++) {
        logd = ff.fastlog( d );
        truelogd = log( d );
        maxdiff = std::max( maxdiff, std::abs( logd-truelogd ) );
        d += 1e-7;
      }
      CPPUNIT_ASSERT( maxdiff < 1.5e-7 );
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

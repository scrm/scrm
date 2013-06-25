#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../src/param.h"

class TestParam : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestParam );

  CPPUNIT_TEST( testParse );
  
  CPPUNIT_TEST_SUITE_END();


 public:
  void testParse() {
    char *argv[] = { "A", "BCD", "EFG" };
    Param pars = Param();
    pars.init();
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestParam );

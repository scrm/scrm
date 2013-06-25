#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/lexical_cast.hpp> 

#include "../src/param.h"
#include "../src/model.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestParam : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestParam );

  CPPUNIT_TEST( testParse );
  CPPUNIT_TEST( testReadInput );
  
  CPPUNIT_TEST_SUITE_END();


 public:
  void testParse() {
    CPPUNIT_ASSERT_NO_THROW( Param().parse() );

    char *argv[] = { "scrm", "15", "10", "-t", "3.74", "-r", "1.24", "1024", "-I", "3", "7", "8", "5" };
    Model model = Param(13, argv).parse();
    CPPUNIT_ASSERT_EQUAL( (size_t)15, model.total_sample_size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)10, model.loci_number() );
    CPPUNIT_ASSERT_EQUAL( (double)3.74, model.mutation_rate() );
    CPPUNIT_ASSERT_EQUAL( (double)1.24, model.recombination_rate() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1024, model.loci_length() );
    CPPUNIT_ASSERT_EQUAL( (size_t)3, model.population_number() );
    CPPUNIT_ASSERT_EQUAL( (size_t)7, model.sample_size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)8, model.sample_size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)5, model.sample_size(2) );
  }

  void testReadInput() {
    CPPUNIT_ASSERT_EQUAL( (int)1, readInput<int>("1") );
    CPPUNIT_ASSERT_EQUAL( (size_t)7, readInput<size_t>("7") );
    CPPUNIT_ASSERT_EQUAL( (double)3.1, readInput<double>("3.1") );
    CPPUNIT_ASSERT_THROW( readInput<int>("ABC"), boost::bad_lexical_cast );
    CPPUNIT_ASSERT_THROW( readInput<int>("-I"), boost::bad_lexical_cast );
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestParam );

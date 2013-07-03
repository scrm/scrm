#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <boost/lexical_cast.hpp> 

#include "../src/param.h"
#include "../src/model.h"
#include "../src/forest.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestParam : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestParam );

  CPPUNIT_TEST( testParse );
  CPPUNIT_TEST( testReadInput );
  
  CPPUNIT_TEST_SUITE_END();

 public:
  void testParse() {
    CPPUNIT_ASSERT_NO_THROW( Param().parse() );

    char *argv[] = { "scrm", "4", "7", "-t", "4.004", "-r", "1.24", "1001", "-l", "1000", "-seed", "123"};
    Param pars = Param(12, argv);
    Model* model = pars.parse();
    //CPPUNIT_ASSERT_EQUAL( (size_t)4, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)7, model->loci_number() );
    CPPUNIT_ASSERT_EQUAL( (double)1e-07, model->mutation_rate() );
    CPPUNIT_ASSERT_EQUAL( (double)3.1e-08, model->recombination_rate() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1001, model->loci_length() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1000, model->exact_window_length() );
    CPPUNIT_ASSERT_EQUAL( (int)123, pars.random_seed );
    CPPUNIT_ASSERT_EQUAL( false, pars.tree_bool );
    delete model;

    char *argv2[] = { "scrm", "15", "10", "-t", "3.74", "-I", "3", "7", "8", "5" };
    CPPUNIT_ASSERT_THROW( Param(10, argv2).parse(), std::invalid_argument ); 
    argv2[3] = "-tv";
    CPPUNIT_ASSERT_THROW( Param(10, argv2).parse(), std::invalid_argument ); 

    char *argv3[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "-T" };
    Param pars2 = Param(11, argv3);
    model = pars2.parse(); 
    CPPUNIT_ASSERT_EQUAL( (size_t)3, model->population_number() );
    CPPUNIT_ASSERT_EQUAL( model->sample_population(4), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model->sample_population(10), (size_t)1 );
    CPPUNIT_ASSERT_EQUAL( model->sample_population(17), (size_t)2 );
    CPPUNIT_ASSERT_EQUAL( model->sample_time(4), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( model->sample_time(17), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( true, pars2.tree_bool );
    delete model;

    char *argv4[] = { "scrm", "23", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "-eI", "12.3", "2", "0", "1" };
    CPPUNIT_ASSERT_NO_THROW( model = Param(15, argv4).parse() ); 
    CPPUNIT_ASSERT_EQUAL( model->sample_population(20), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model->sample_population(22), (size_t)2 );
    CPPUNIT_ASSERT_EQUAL( model->sample_time(20), (double)12.3 );
    CPPUNIT_ASSERT_EQUAL( model->sample_time(22), (double)12.3 );
    delete model;
    
    char *argv5[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "-eN", "8.2", "0.75" };
    model = Param(13, argv5).parse(); 
    model->resetTime();
    model->increaseTime();
    CPPUNIT_ASSERT_EQUAL( model->population_size(0), (size_t)(0.75*model->default_pop_size) );
    CPPUNIT_ASSERT_EQUAL( model->population_size(1), (size_t)(0.75*model->default_pop_size) );
    CPPUNIT_ASSERT_EQUAL( model->population_size(2), (size_t)(0.75*model->default_pop_size) );
    delete model;
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

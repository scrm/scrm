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
  CPPUNIT_TEST( testParseMigrationOptions );
  CPPUNIT_TEST( testParseGrowthOptions );
  CPPUNIT_TEST( testReadInput );
  
  CPPUNIT_TEST_SUITE_END();

 private:
    Model* model;

 public:
  void testParse() {
    CPPUNIT_ASSERT_NO_THROW( Param().parse() );

    char *argv[] = { "scrm", "4", "7", "-t", "4.004", "-r", "1.24", "1001", "-l", "1000", "-seed", "123"};
    Param pars = Param(12, argv);
    model = pars.parse();
    //CPPUNIT_ASSERT_EQUAL( (size_t)4, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)7, model->loci_number() );
    CPPUNIT_ASSERT_EQUAL( 4.004/(4*model->default_pop_size), model->mutation_rate() );
    CPPUNIT_ASSERT_EQUAL( 1.24/(4*model->default_pop_size*1000), model->rec_rate() );
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
   
    // -N & -eN 
    char *argv5[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", 
                      "-N", "0.3", "-eN", "8.2", "0.75", "-G", "1.5" };
    CPPUNIT_ASSERT_NO_THROW( model = Param(17, argv5).parse() ); 
    model->resetTime();
    CPPUNIT_ASSERT_EQUAL( 0.3*model->default_pop_size, model->population_size(0) );
    CPPUNIT_ASSERT_EQUAL( 0.3*model->default_pop_size, model->population_size(1) );
    CPPUNIT_ASSERT_EQUAL( 0.3*model->default_pop_size, model->population_size(2) );
    model->increaseTime();
    CPPUNIT_ASSERT_EQUAL( 0.75*model->default_pop_size, model->population_size(0) );
    CPPUNIT_ASSERT_EQUAL( 0.75*model->default_pop_size, model->population_size(1) );
    CPPUNIT_ASSERT_EQUAL( 0.75*model->default_pop_size, model->population_size(2) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model->growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model->growth_rate(1) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model->growth_rate(2) );
    delete model;

    // -n & -en 
    char *argv6[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "-G", "1.5", 
                      "-n", "1", "0.3", "-eN", "1.1", "0.75", "-en", "2", "2", "0.1", "-eG", "1.5", "2" };
    CPPUNIT_ASSERT_NO_THROW( model = Param(25, argv6).parse() ); 
    model->finalize();
    model->resetTime();
    CPPUNIT_ASSERT_EQUAL( (double)model->default_pop_size, model->population_size(0) );
    CPPUNIT_ASSERT_EQUAL( 0.3*model->default_pop_size, model->population_size(1) );
    CPPUNIT_ASSERT_EQUAL( (double)model->default_pop_size, model->population_size(2) );
    CPPUNIT_ASSERT_EQUAL( 1.5, model->growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model->growth_rate(1) );
    CPPUNIT_ASSERT_EQUAL( 1.5, model->growth_rate(2) );
    model->increaseTime();
    model->increaseTime();
    model->increaseTime();
    CPPUNIT_ASSERT_EQUAL( 0.75*model->default_pop_size, model->population_size(0) );
    CPPUNIT_ASSERT_EQUAL( 0.75*model->default_pop_size, model->population_size(1) );
    CPPUNIT_ASSERT_EQUAL( 0.10*model->default_pop_size, model->population_size(2) );
    CPPUNIT_ASSERT_EQUAL( 2.0, model->growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 2.0, model->growth_rate(1) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model->growth_rate(2) );
    delete model;
  }

  void testParseMigrationOptions() {
    // -ma
    char *argv[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
                      "-ma", "x", "5", "7", "x" };
    CPPUNIT_ASSERT_NO_THROW( model = Param(14, argv).parse(); );
    model->resetTime();
    CPPUNIT_ASSERT_EQUAL( 5.0/(4*model->default_pop_size), model->migration_rate(0, 1) );
    CPPUNIT_ASSERT_EQUAL( 7.0/(4*model->default_pop_size), model->migration_rate(1, 0) );
    delete model;

    // -ema
    char *argv2[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
                      "-ema", "1.6", "x", "5", "7", "x" };
    CPPUNIT_ASSERT_NO_THROW( model = Param(15, argv2).parse(); );
    model->resetTime();
    model->increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.6, model->getCurrentTime() );
    CPPUNIT_ASSERT_EQUAL( 5.0/(4*model->default_pop_size), model->migration_rate(0, 1) );
    CPPUNIT_ASSERT_EQUAL( 7.0/(4*model->default_pop_size), model->migration_rate(1, 0) );
    delete model;

    // -M
    char *argv3[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "10", "10", "0", 
                      "-M", "5" };
    CPPUNIT_ASSERT_NO_THROW( model = Param(12, argv3).parse(); );
    model->resetTime();
    CPPUNIT_ASSERT_EQUAL( 2.5/(4*model->default_pop_size), model->migration_rate(1, 0) );
    CPPUNIT_ASSERT_EQUAL( 2.5/(4*model->default_pop_size), model->migration_rate(0, 1) );
    CPPUNIT_ASSERT_EQUAL( 2.5/(4*model->default_pop_size), model->migration_rate(0, 2) );
    delete model;

    // -eM
    char *argv4[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "10", "10", "0", 
                      "-eM", "1.6", "5" };
    CPPUNIT_ASSERT_NO_THROW( model = Param(13, argv4).parse(); );
    model->resetTime();
    model->increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.6, model->getCurrentTime() );
    CPPUNIT_ASSERT_EQUAL( 2.5/(4*model->default_pop_size), model->migration_rate(1, 0) );
    CPPUNIT_ASSERT_EQUAL( 2.5/(4*model->default_pop_size), model->migration_rate(0, 1) );
    CPPUNIT_ASSERT_EQUAL( 2.5/(4*model->default_pop_size), model->migration_rate(0, 2) );
    delete model;

    // -esme
    char *argv5[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
                      "-esme", "1.6", "1", "0", "0.5", "-esme", "1.6", "0", "1", "0.1" };
    CPPUNIT_ASSERT_NO_THROW( model = Param(19, argv5).parse(); );
    model->resetTime();
    model->increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.6, model->getCurrentTime() );
    CPPUNIT_ASSERT_EQUAL( 0.1, model->single_mig_pop(0, 1) );
    CPPUNIT_ASSERT_EQUAL( 0.5, model->single_mig_pop(1, 0) );
    delete model;

    // -ej
    char *argv6[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
                      "-M", "1.3", "-ej", "1.6", "1", "0" };
    CPPUNIT_ASSERT_NO_THROW( model = Param(15, argv6).parse(); );
    model->resetTime();
    model->increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.6, model->getCurrentTime() );
    CPPUNIT_ASSERT_EQUAL( 1.0, model->single_mig_pop(1, 0) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model->single_mig_pop(0, 1) );

    CPPUNIT_ASSERT_EQUAL( 0.0, model->migration_rate(0, 1) );
    CPPUNIT_ASSERT_EQUAL( 1.3/(4*model->default_pop_size), model->migration_rate(1, 0) );
    delete model;
  }

  void testReadInput() {
    CPPUNIT_ASSERT_EQUAL( (int)1, readInput<int>("1") );
    CPPUNIT_ASSERT_EQUAL( (size_t)7, readInput<size_t>("7") );
    CPPUNIT_ASSERT_EQUAL( (double)3.1, readInput<double>("3.1") );
    CPPUNIT_ASSERT_THROW( readInput<int>("ABC"), boost::bad_lexical_cast );
    CPPUNIT_ASSERT_THROW( readInput<int>("-I"), boost::bad_lexical_cast );
  }

  void testParseGrowthOptions() {
    // -G && -eG
    char *argv[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
      "-G", "2", "-eG", "1", "3" };
    CPPUNIT_ASSERT_NO_THROW( model = Param(14, argv).parse(); );
    model->resetTime();
    CPPUNIT_ASSERT_EQUAL( 2.0, model->growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 2.0, model->growth_rate(1) );
    model->increaseTime();
    CPPUNIT_ASSERT_EQUAL( 3.0, model->growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 3.0, model->growth_rate(1) );
    delete model;

    // -g && -eg
    char *argv2[] = { 
      "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
      "-g", "1", "0.1", "-eG", "1", "3", "-eg", "2", "0", "2.4" };
    CPPUNIT_ASSERT_NO_THROW( model = Param(19, argv2).parse(); );
    model->resetTime();
    CPPUNIT_ASSERT_EQUAL( model->default_growth_rate, model->growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 0.1, model->growth_rate(1) );
    model->increaseTime();
    model->increaseTime();
    CPPUNIT_ASSERT_EQUAL( 2.4, model->growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 3.0, model->growth_rate(1) );
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestParam );

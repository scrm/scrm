#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../../src/param.h"
#include "../../src/model.h"
#include "../../src/forest.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestParam : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestParam );

  CPPUNIT_TEST( testStringConstructor );
  CPPUNIT_TEST( testParseBasicCmdStructure );
  CPPUNIT_TEST( testParseSeeds );
  CPPUNIT_TEST( testParseCommonOptions );
  CPPUNIT_TEST( testParseMigrationOptions );
  CPPUNIT_TEST( testParseOutputOptions );
  CPPUNIT_TEST( testParseGrowthOptions );
  CPPUNIT_TEST( testParseVariableRates );
  CPPUNIT_TEST( testParseUnsupportedMsArguments );
  CPPUNIT_TEST( testParseSplitOptions );
  CPPUNIT_TEST( testParseMergeOptions );
  CPPUNIT_TEST( testParseSequenceScaling );
  CPPUNIT_TEST( testErrorOnInfintieRecRate );
  CPPUNIT_TEST( testScientificNotation );

  CPPUNIT_TEST_SUITE_END();

 private:
  Model model;

 public:
  void testParseBasicCmdStructure() {
    Param pars;
    char *argv1[] = { "scrm" };
    CPPUNIT_ASSERT_THROW( Param(1, argv1).parse(model), std::invalid_argument ); 

    char *argv2[] = { "scrm", "-h" };
    pars = Param(2, argv2);
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT( pars.help() ); 

    char *argv3[] = { "scrm", "--help" };
    pars = Param(2, argv3);
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT( pars.help() ); 

    char *argv4[] = { "scrm", "2", "1", "-h" };
    pars = Param(4, argv4);
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT( pars.help() ); 

    char *argv5[] = { "scrm", "1" };
    CPPUNIT_ASSERT_THROW( Param(2, argv5).parse(model), std::invalid_argument ); 

    char *argv6[] = { "scrm", "1", "-t", "5" };
    CPPUNIT_ASSERT_THROW( Param(4, argv6).parse(model), std::invalid_argument ); 

    char *argv7[] = { "scrm", "-t", "5" };
    CPPUNIT_ASSERT_THROW( Param(3, argv7).parse(model), std::invalid_argument ); 

    char *argv8[] = { "scrm", "2", "1", "-t", "5" };
    pars = Param(5, argv8);
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT( !pars.help() ); 
    CPPUNIT_ASSERT( !pars.version() ); 

    char *argv9[] = { "scrm", "-v" };
    pars = Param(2, argv9);
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT( pars.version() ); 

    char *argv10[] = { "scrm", "--version" };
    pars = Param(2, argv10);
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT( pars.version() ); 

    char *argv11[] = { "scrm", "2", "1", "-v" };
    pars = Param(4, argv11);
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT( pars.version() ); 
  }

  void testParseCommonOptions() {
    CPPUNIT_ASSERT_NO_THROW( Param().parse(model) );

    char *argv[] = { "scrm", "4", "7", "-t", "40.04", 
      "-r", "1.24", "1001", "-l", "1000"};

    Param pars = Param(10, argv);
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT( !pars.help() );

    CPPUNIT_ASSERT_EQUAL( (size_t)4, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)7, model.loci_number() );
    CPPUNIT_ASSERT( areSame(40.04/(4*model.default_pop_size*1001), model.mutation_rate()) );
    CPPUNIT_ASSERT( areSame(1.24/(4*model.default_pop_size*1000), model.recombination_rate()) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1001, model.loci_length() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1000, model.exact_window_length() );

    char *argv2[] = { "scrm", "15", "10", "-t", "3.74", "-I", "3", "7", "8", "5" };
    CPPUNIT_ASSERT_THROW( Param(10, argv2).parse(model), std::invalid_argument ); 
    argv2[3] = "-tv";
    CPPUNIT_ASSERT_THROW( Param(10, argv2).parse(model), std::invalid_argument ); 

    char *argv3[] = { "scrm", "20", "10", "-t", "3.74", 
                      "-I", "3", "7", "8", "5", "-T", "-M", "5.0" };
    Param pars2 = Param(13, argv3);
    CPPUNIT_ASSERT_NO_THROW( pars2.parse(model) ); 
    CPPUNIT_ASSERT( !pars2.help() );
    CPPUNIT_ASSERT_EQUAL( (size_t)3, model.population_number() );
    CPPUNIT_ASSERT_EQUAL( (size_t)20, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(4), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(10), (size_t)1 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(17), (size_t)2 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(4), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(17), (double)0.0 );

    char *argv32[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "5.0", "-T"};
    CPPUNIT_ASSERT_NO_THROW( Param(12, argv32).parse(model) ); 
    CPPUNIT_ASSERT_EQUAL( (size_t)3, model.population_number() );
    CPPUNIT_ASSERT_EQUAL( (size_t)20, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(4), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(10), (size_t)1 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(17), (size_t)2 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(4), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(17), (double)0.0 );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 2)) );

    char *argv33[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "5.0"};
    CPPUNIT_ASSERT_NO_THROW( Param(11, argv33).parse(model) ); 
    CPPUNIT_ASSERT_EQUAL( (size_t)3, model.population_number() );
    CPPUNIT_ASSERT_EQUAL( (size_t)20, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(4), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(10), (size_t)1 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(17), (size_t)2 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(4), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(17), (double)0.0 );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 2)) );

    char *argv34[] = { "scrm", "2", "1", "-t", "3.74", "-I", "2", "1", "1", "5.0"};
    CPPUNIT_ASSERT_NO_THROW( Param(10, argv34).parse(model) ); 
    CPPUNIT_ASSERT_EQUAL( (size_t)2, model.population_number() );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(0), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(1), (size_t)1 );

    char *argv4[] = { "scrm", "23", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "-eI", "12.3", "2", "0", "1" , "-M", "5.0" };
    CPPUNIT_ASSERT_NO_THROW( Param(17, argv4).parse(model) ); 
    CPPUNIT_ASSERT_EQUAL( model.sample_population(20), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(22), (size_t)2 );
    CPPUNIT_ASSERT_EQUAL( 12.3 * 4 * model.default_pop_size, model.sample_time(20) );
    CPPUNIT_ASSERT_EQUAL( 12.3 * 4 * model.default_pop_size, model.sample_time(22) );

    // -N & -eN 
    char *argv5[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", 
      "-N", "0.3", "-eN", "8.2", "0.75", "-G", "1.5", "-M", "5.0"};
    CPPUNIT_ASSERT_NO_THROW( Param(19, argv5).parse(model) ); 
    model.resetTime();
    CPPUNIT_ASSERT_EQUAL( 0.0, model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame(0.3*model.default_pop_size, model.population_size(0)) );
    CPPUNIT_ASSERT( areSame(0.3*model.default_pop_size, model.population_size(1)) );
    CPPUNIT_ASSERT( areSame(0.3*model.default_pop_size, model.population_size(2)) );
    model.increaseTime();
    CPPUNIT_ASSERT( areSame(8.2 * 4 * model.default_pop_size, model.getCurrentTime()) );
    CPPUNIT_ASSERT( areSame(0.75*model.default_pop_size, model.population_size(0)) );
    CPPUNIT_ASSERT( areSame(0.75*model.default_pop_size, model.population_size(1)) );
    CPPUNIT_ASSERT( areSame(0.75*model.default_pop_size, model.population_size(2)) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate(1) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate(2) );

    // -n & -en 
    char *argv6[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "7", "8", "5", "-G", "1.5", 
      "-n", "2", "0.3", "-eN", "1.1", "0.75", "-en", "2", "3", "0.1", "-eG", "1.5", "2", "-M", "5.0" };
    CPPUNIT_ASSERT_NO_THROW( Param(27, argv6).parse(model) ); 
    model.finalize();
    //std::cout << model << std::endl;
    model.resetTime();
    CPPUNIT_ASSERT_EQUAL( 0.0, model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame(model.default_pop_size, model.population_size(0)) );
    CPPUNIT_ASSERT( areSame(0.3*model.default_pop_size, model.population_size(1)) );
    CPPUNIT_ASSERT( areSame(model.default_pop_size, (size_t)model.population_size(2)) );
    CPPUNIT_ASSERT( areSame(1.5 / 4 / model.default_pop_size, model.growth_rate(0)) );
    CPPUNIT_ASSERT( areSame(1.5 / 4 / model.default_pop_size, model.growth_rate(1)) );
    CPPUNIT_ASSERT( areSame(1.5 / 4 / model.default_pop_size, model.growth_rate(2)) );
    model.increaseTime();
    CPPUNIT_ASSERT( areSame(1.1 * 4 * model.default_pop_size, model.getCurrentTime()) );
    model.increaseTime();
    CPPUNIT_ASSERT( areSame(1.5 * 4 * model.default_pop_size, model.getCurrentTime()) );
    model.increaseTime();
    CPPUNIT_ASSERT( areSame(2.0 * 4 * model.default_pop_size, model.getCurrentTime()) );
    CPPUNIT_ASSERT( 0.75*model.default_pop_size > model.population_size(0) );
    CPPUNIT_ASSERT( 0.75*model.default_pop_size > model.population_size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)(0.10*model.default_pop_size), (size_t)model.population_size(2) );
    CPPUNIT_ASSERT_EQUAL( 2.0 / 4 / model.default_pop_size  , model.growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 2.0 / 4 / model.default_pop_size, model.growth_rate(1) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate(2) );
  }

  void testParseMigrationOptions() {
    // -ma
    char *argv[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
      "-ma", "x", "5", "7", "x" };
    CPPUNIT_ASSERT_NO_THROW( Param(14, argv).parse(model); );
    model.resetTime();
    CPPUNIT_ASSERT( areSame(5.0/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(7.0/(4*model.default_pop_size), model.migration_rate(1, 0)) );

    // -ema
    char *argv2[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
      "-ema", "1.6", "x", "5", "7", "x" };
    CPPUNIT_ASSERT_NO_THROW( Param(15, argv2).parse(model); );
    model.resetTime();
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.6 * 4 * model.default_pop_size, model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame(5.0/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(7.0/(4*model.default_pop_size), model.migration_rate(1, 0)) );

    // -M
    char *argv3[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "10", "10", "0", 
      "-M", "5" };
    CPPUNIT_ASSERT_NO_THROW( Param(12, argv3).parse(model); );
    model.resetTime();
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 2)) );

    // -eM
    char *argv4[] = { "scrm", "20", "10", "-t", "3.74", "-I", "3", "10", "10", "0", 
      "-eM", "1.6", "5" };
    CPPUNIT_ASSERT_NO_THROW( Param(13, argv4).parse(model); );
    model.resetTime();
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.6 * 4 * model.default_pop_size,  model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 2)) );

    CPPUNIT_ASSERT_NO_THROW( Param("scrm 20 1 -I 2 13 7 -m 1 2 1.5").parse(model) );
    CPPUNIT_ASSERT( areSame(1.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(0.0/(4*model.default_pop_size), model.migration_rate(1, 0)) );

    CPPUNIT_ASSERT_NO_THROW( Param("scrm 20 1 -I 2 13 7 -m 1 2 1.5 -m 2 1 0.5").parse(model) );
    CPPUNIT_ASSERT( areSame(1.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(0.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );

    CPPUNIT_ASSERT_NO_THROW( Param("scrm 20 1 -I 2 13 7 -em 1.0 2 1 0.5 -em 2.0 2 1 0.6").parse(model) );
    CPPUNIT_ASSERT( areSame(0.0/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    CPPUNIT_ASSERT( areSame(0.0/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    model.increaseTime();
    CPPUNIT_ASSERT( areSame(0.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    CPPUNIT_ASSERT( areSame(0.0/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    model.increaseTime();
    CPPUNIT_ASSERT( areSame(0.6/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    CPPUNIT_ASSERT( areSame(0.0/(4*model.default_pop_size), model.migration_rate(0, 1)) );
  }

  void testParseMergeOptions() {
    // -es
    char *argv1[] = { "scrm", "20", "10", 
                      "-I", "2", "10", "10", "1.5", 
                      "-es", "1.6", "2", "0.5" };
    CPPUNIT_ASSERT_NO_THROW( Param(12, argv1).parse(model); );
    model.resetTime();
    CPPUNIT_ASSERT_EQUAL( (size_t)3, model.population_number() );
    CPPUNIT_ASSERT( areSame(1.5/(4*model.default_pop_size), model.migration_rate(0,1)) );
    CPPUNIT_ASSERT( areSame(1.5/(4*model.default_pop_size), model.migration_rate(1,0)) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.migration_rate(0,2) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.migration_rate(2,0) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.migration_rate(0,2) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.migration_rate(2,1) );
    CPPUNIT_ASSERT_EQUAL( model.default_pop_size, model.population_size(2) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate(2) );

    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.6 * 4 * model.default_pop_size, model.getCurrentTime() );
    CPPUNIT_ASSERT_EQUAL( 0.5, model.single_mig_pop(1, 2) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.single_mig_pop(2, 1) );
    CPPUNIT_ASSERT( areSame(1.5/(4*model.default_pop_size), model.migration_rate(0,1)) );
    CPPUNIT_ASSERT( areSame(1.5/(4*model.default_pop_size), model.migration_rate(1,0)) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.migration_rate(0,2) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.migration_rate(2,0) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.migration_rate(1,2) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.migration_rate(2,1) );
    CPPUNIT_ASSERT_EQUAL( model.default_pop_size, model.population_size(1) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate() );

    char *argv2[] = { "scrm", "20", "10", 
                      "-I", "2", "10", "10", "1.5", 
                      "-es", "1.6", "2", "0.5",
                      "-eM", "0.9", "5.0" };
    CPPUNIT_ASSERT_THROW( Param(15, argv2).parse(model), std::invalid_argument );

    char *argv3[] = { "scrm", "20", "10", 
                      "-I", "2", "10", "10", "1.5", 
                      "-es", "1.6", "2", "0.5",
                      "-eG", "0.9", "5.0" };
    CPPUNIT_ASSERT_THROW( Param(15, argv3).parse(model), std::invalid_argument );

    char *argv4[] = { "scrm", "20", "10", 
                      "-I", "2", "10", "10", "1.5", 
                      "-es", "1.6", "2", "0.5",
                      "-eN", "0.9", "5.0" };
    CPPUNIT_ASSERT_THROW( Param(15, argv4).parse(model), std::invalid_argument );

    char *argv5[] = { "scrm", "20", "10", 
                      "-I", "2", "10", "10", "1.5", 
                      "-es", "1.6", "2", "0.5",
                      "-es", "0.9", "3", "1" };
    CPPUNIT_ASSERT_THROW( Param(16, argv5).parse(model), std::invalid_argument );
  }

  void testParseSplitOptions() {
    // -ej
    char *argv6[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
      "-ej", "1.6", "2", "1", "-M", "1.3" };
    CPPUNIT_ASSERT_NO_THROW( Param(15, argv6).parse(model); );
    model.resetTime();
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.6 * 4 * model.default_pop_size, model.getCurrentTime() );
    CPPUNIT_ASSERT_EQUAL( 1.0, model.single_mig_pop(1, 0) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.single_mig_pop(0, 1) );

    CPPUNIT_ASSERT_EQUAL( 0.0, model.migration_rate(0, 1) );
    CPPUNIT_ASSERT_EQUAL( 1.3/(4*model.default_pop_size), model.migration_rate(1, 0) );
  }

  void testParseOutputOptions() {
    char *argv[] = { "scrm", "20", "10", "-T" };
    CPPUNIT_ASSERT_NO_THROW( Param(4, argv).parse(model); );
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 1 );

    char *argv2[] = { "scrm", "20", "10", "-O" };
    CPPUNIT_ASSERT_NO_THROW( Param(4, argv2).parse(model); );
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 1 );

    char *argv3[] = { "scrm", "20", "10", "-T", "-O" };
    CPPUNIT_ASSERT_NO_THROW( Param(5, argv3).parse(model); );
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 2 );

    char *argv4[] = { "scrm", "20", "10", "-t", "5" };
    CPPUNIT_ASSERT_NO_THROW( Param(5, argv4).parse(model); );
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 1 );

    char *argv5[] = { "scrm", "20", "10", "-L" };
    CPPUNIT_ASSERT_NO_THROW( Param(4, argv5).parse(model); );
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 1 );

    char *argv6[] = { "scrm", "20", "10", "-oSFS" };
    CPPUNIT_ASSERT_THROW( Param(4, argv6).parse(model), std::invalid_argument ); 

    char *argv7[] = { "scrm", "20", "10", "-oSFS", "-t", "5" };
    CPPUNIT_ASSERT_NO_THROW( Param(6, argv7).parse(model); );
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 2 );

    char *argv8[] = { "scrm", "20", "10", "-oSFS", "-t", "5", "-L" };
    CPPUNIT_ASSERT_NO_THROW( Param(7, argv8).parse(model); );
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 3 );

    char *argv9[] = { "scrm", "20", "10", "-oSFS", "-t", "5", "-L", "-T" };
    CPPUNIT_ASSERT_NO_THROW( Param(8, argv9).parse(model); );
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 4 );
  }

  void testParseGrowthOptions() {
    // -G && -eG
    char *argv[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
      "-G", "2", "-eG", "1", "3", "-M", "5.0" };
    CPPUNIT_ASSERT_NO_THROW( Param(16, argv).parse(model); );
    model.resetTime();
    CPPUNIT_ASSERT( areSame(2.0 / 4 / model.default_pop_size, model.growth_rate(0)) );
    CPPUNIT_ASSERT( areSame(2.0 / 4 / model.default_pop_size, model.growth_rate(1)) );
    model.increaseTime();
    CPPUNIT_ASSERT( areSame(1.0 * 4 * model.default_pop_size, model.getCurrentTime()) );
    CPPUNIT_ASSERT( areSame(3.0 / 4 / model.default_pop_size, model.growth_rate(0)) );
    CPPUNIT_ASSERT( areSame(3.0 / 4 / model.default_pop_size, model.growth_rate(1)) );

    // -g && -eg
    char *argv2[] = { 
      "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
      "-g", "2", "0.1", "-eG", "1", "3", "-eg", "2", "1", "2.4", "-M", "5.0"};
    CPPUNIT_ASSERT_NO_THROW( Param(21, argv2).parse(model); );
    model.resetTime();
    CPPUNIT_ASSERT_EQUAL( model.default_growth_rate, model.growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 0.1 / 4 / model.default_pop_size, model.growth_rate(1) );
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.0 * 4 * model.default_pop_size, model.getCurrentTime() );
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 2.0 * 4 * model.default_pop_size, model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame(2.4 / 4 / model.default_pop_size, model.growth_rate(0)) );
    CPPUNIT_ASSERT( areSame(3.0 / 4 / model.default_pop_size, model.growth_rate(1)) );
  }

  void testParseVariableRates() {
    // -st
    char *argv[] = { "scrm", "20", "10", "-t", "3.74", "-st", "10", "5.1", "-st", "4.5", "3.2"};
    CPPUNIT_ASSERT_NO_THROW( Param(11, argv).parse(model); );
    double scale = 1 / ( 4 * model.default_pop_size * model.default_loci_length );
    CPPUNIT_ASSERT( areSame(3.74 * scale, model.mutation_rate()) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.getCurrentSequencePosition() );
    CPPUNIT_ASSERT_EQUAL( 4.5, model.getNextSequencePosition() );

    model.increaseSequencePosition();
    CPPUNIT_ASSERT( areSame(3.2 * scale, model.mutation_rate()) );
    CPPUNIT_ASSERT_EQUAL( 4.5, model.getCurrentSequencePosition() );
    CPPUNIT_ASSERT_EQUAL( 10.0, model.getNextSequencePosition() );

    model.increaseSequencePosition();
    CPPUNIT_ASSERT( areSame(5.1 * scale, model.mutation_rate()) );
    CPPUNIT_ASSERT_EQUAL( 10.0, model.getCurrentSequencePosition() );

    // -sr
    char *argv2[] = { "scrm", "20", "10", "-r", "3.74", "100", "-sr", "10", "5.1", "-sr", "4.5", "3.2"};
    CPPUNIT_ASSERT_NO_THROW( Param(12, argv2).parse(model); );
    model.resetSequencePosition();
    scale = 1 / ( 4 * model.default_pop_size * 99 );
    CPPUNIT_ASSERT( areSame(3.74 * scale, model.recombination_rate()) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.getCurrentSequencePosition() );
    CPPUNIT_ASSERT_EQUAL( 4.5, model.getNextSequencePosition() );

    model.increaseSequencePosition();
    CPPUNIT_ASSERT( areSame(3.2 * scale, model.recombination_rate()) );
    CPPUNIT_ASSERT_EQUAL( 4.5, model.getCurrentSequencePosition() );
    CPPUNIT_ASSERT_EQUAL( 10.0, model.getNextSequencePosition() );

    model.increaseSequencePosition();
    CPPUNIT_ASSERT( areSame(5.1 * scale, model.recombination_rate()) );
    CPPUNIT_ASSERT_EQUAL( 10.0, model.getCurrentSequencePosition() );

    // both
    char *argv3[] = { "scrm", "20", "10", "-r", "3.74", "100", "-sr", "10", "5.1", "-st", "4.5", "3.2"};
    CPPUNIT_ASSERT_NO_THROW( Param(12, argv3).parse(model); );
    model.resetSequencePosition();
    model.resetTime();
    scale = 1 / ( 4 * model.default_pop_size * 99 );
    double scale2 = 1 / ( 4 * model.default_pop_size * 100 );

    CPPUNIT_ASSERT( areSame(3.74 * scale, model.recombination_rate()) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.mutation_rate() );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.getCurrentSequencePosition() );
    CPPUNIT_ASSERT_EQUAL( 4.5, model.getNextSequencePosition() );

    model.increaseSequencePosition();
    CPPUNIT_ASSERT_EQUAL( 4.5, model.getCurrentSequencePosition() );
    CPPUNIT_ASSERT_EQUAL( 10.0, model.getNextSequencePosition() );
    CPPUNIT_ASSERT( areSame(3.74 * scale, model.recombination_rate()) );
    CPPUNIT_ASSERT( areSame(3.2 * scale2, model.mutation_rate()) );

    model.increaseSequencePosition();
    CPPUNIT_ASSERT( areSame(5.1 * scale, model.recombination_rate()) );
    CPPUNIT_ASSERT( areSame(3.2 * scale2, model.mutation_rate()) );
    CPPUNIT_ASSERT_EQUAL( 10.0, model.getCurrentSequencePosition() );
  }

  void testParseUnsupportedMsArguments() {
    char *argv[] = { "scrm", "20", "10", "-t", "3.74", "-c", "10", "5.1", "-st", "4.5", "3.2"};
    CPPUNIT_ASSERT_THROW( Param(11, argv).parse(model), std::invalid_argument ); 

    char *argv2[] = { "scrm", "20", "10", "-s", "5"};
    CPPUNIT_ASSERT_THROW( Param(5, argv2).parse(model), std::invalid_argument ); 
  }

  void testParseSeeds() {
    Param pars = Param("scrm 4 7 -seed 123 -t 5");
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT_EQUAL( (size_t)123, pars.random_seed() );

    pars = Param("scrm 4 7 -seed 1 2 3 -t 5");
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    size_t seed = pars.random_seed();
    pars = Param("scrm 4 7 -seed 1 2 3 -t 5");
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT_EQUAL( seed, pars.random_seed() );

    char *argv3[] = { "scrm", "20", "10", "-seed", "1", "2", "3", "4"};
    CPPUNIT_ASSERT_THROW( Param(8, argv3).parse(model), std::invalid_argument ); 

    char *argv4[] = { "scrm", "20", "10", "-seed", "-t", "2"};
    CPPUNIT_ASSERT_THROW( Param(6, argv4).parse(model), std::invalid_argument ); 
  }

  void testStringConstructor() {
    Param pars = Param("scrm 4 7 -t 40.04 -r 0.5 1001 -l 1000");
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );

    CPPUNIT_ASSERT_EQUAL( (size_t)4, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)7, model.loci_number() );
    CPPUNIT_ASSERT( areSame(40.04/(4*model.default_pop_size*1001), model.mutation_rate()) );
    CPPUNIT_ASSERT( areSame(0.5/(4*model.default_pop_size*1000), model.recombination_rate()) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1001, model.loci_length() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1000, model.exact_window_length() );
  }

  void testParseSequenceScaling() {
    char *argv1[] = { "scrm", "4", "7", "-SC", "ms"};
    Param pars = Param(5, argv1);
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT( model.getSequenceScaling() == ms );

    char *argv2[] = { "scrm", "4", "7", "-SC", "rel"};
    pars = Param(5, argv2);
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT( model.getSequenceScaling() == relative );

    char *argv3[] = { "scrm", "4", "7", "-SC", "abs"};
    pars = Param(5, argv3);
    CPPUNIT_ASSERT_NO_THROW( pars.parse(model) );
    CPPUNIT_ASSERT( model.getSequenceScaling() == absolute );

    char *argv4[] = { "scrm", "4", "7", "-SC", "blub"};
    CPPUNIT_ASSERT_THROW( Param(5, argv4).parse(model), std::invalid_argument );
    char *argv5[] = { "scrm", "4", "7", "-SC"};
    CPPUNIT_ASSERT_THROW( Param(4, argv5).parse(model), std::invalid_argument );
  }


  void testErrorOnInfintieRecRate() {
    Param pars = Param("scrm 4 7 -r 0.5 1");
    CPPUNIT_ASSERT_THROW(pars.parse(model), std::invalid_argument);

    pars = Param("scrm 4 7 -r 0.5 0");
    CPPUNIT_ASSERT_THROW(pars.parse(model), std::invalid_argument);
  }


  void testScientificNotation() {
    Param pars;
    pars = Param("scrm 4 7 -r 1e3 1001");
    CPPUNIT_ASSERT_NO_THROW(pars.parse(model));
    CPPUNIT_ASSERT_EQUAL(1.0/(4*model.default_pop_size), model.recombination_rate());

    pars = Param("scrm 4 7 -r 0.5 2e3");
    CPPUNIT_ASSERT_NO_THROW(pars.parse(model));
    CPPUNIT_ASSERT_EQUAL((size_t)2000, model.loci_length());

    pars = Param("scrm 4 1e3 -r 0.5 2");
    CPPUNIT_ASSERT_NO_THROW(pars.parse(model));
    CPPUNIT_ASSERT_EQUAL((size_t)1000, model.loci_number());
  }

};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestParam );

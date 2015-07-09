#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../../src/param.h"
#include "../../src/model.h"
#include "../../src/forest.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestParam : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestParam );

  CPPUNIT_TEST( testMainConstructor );
  CPPUNIT_TEST( testStringConstructor );
  CPPUNIT_TEST( testParseBasicCmdStructure );
  CPPUNIT_TEST( testParseSeeds );
  CPPUNIT_TEST( testParseCommonOptions );
  CPPUNIT_TEST( testParseMigrationOptions );
  CPPUNIT_TEST( testParseMergeOptions );
  CPPUNIT_TEST( testParseSplitOptions );
  CPPUNIT_TEST( testParseOutputOptions );
  CPPUNIT_TEST( testParseGrowthOptions );
  CPPUNIT_TEST( testParseVariableRates );
  CPPUNIT_TEST( testParseUnsupportedMsArguments );
  CPPUNIT_TEST( testParseSequenceScaling );
  CPPUNIT_TEST( testErrorOnInfintieRecRate );
  CPPUNIT_TEST( testScientificNotation );
  CPPUNIT_TEST( testApproximation );

  CPPUNIT_TEST_SUITE_END();

 private:
  Model model;

 public:
  void testMainConstructor() {
    char *argv[] = { "scrm", "arg1", "arg2", "arg3" };
    Param pars(4, argv);
    CPPUNIT_ASSERT_EQUAL((size_t)3, pars.argv_.size());
    CPPUNIT_ASSERT_EQUAL(std::string("arg1"), pars.argv_[0]);
    CPPUNIT_ASSERT_EQUAL(std::string("arg2"), pars.argv_[1]);
    CPPUNIT_ASSERT_EQUAL(std::string("arg3"), pars.argv_[2]);
  }

  void testStringConstructor() {
    Param pars("arg1 arg2 arg3");
    CPPUNIT_ASSERT_EQUAL((size_t)3, pars.argv_.size());
    CPPUNIT_ASSERT_EQUAL(std::string("arg1"), pars.argv_[0]);
    CPPUNIT_ASSERT_EQUAL(std::string("arg2"), pars.argv_[1]);
    CPPUNIT_ASSERT_EQUAL(std::string("arg3"), pars.argv_[2]);
  }

  void testParseBasicCmdStructure() {
    Model model;
    CPPUNIT_ASSERT_THROW(Param("scrm").parse(), std::invalid_argument); 

    Param pars = Param("-h");
    model = pars.parse();
    CPPUNIT_ASSERT( pars.help() ); 

    pars = Param("--help");
    model = pars.parse();
    CPPUNIT_ASSERT( pars.help() ); 

    pars = Param("2 1 -h");
    model = pars.parse();
    CPPUNIT_ASSERT(pars.help()); 

    CPPUNIT_ASSERT_THROW( Param("1").parse(), std::invalid_argument ); 
    CPPUNIT_ASSERT_THROW( Param("1 -t 5").parse(), std::invalid_argument ); 
    CPPUNIT_ASSERT_THROW( Param("-t 5").parse(), std::invalid_argument ); 

    pars = Param("2 1");
    model = pars.parse();
    CPPUNIT_ASSERT( !pars.help() ); 
    CPPUNIT_ASSERT( !pars.version() ); 

    pars = Param("-v");
    model = pars.parse();
    CPPUNIT_ASSERT( pars.version() ); 

    pars = Param("--version");
    model = pars.parse();
    CPPUNIT_ASSERT( pars.version() ); 

    pars = Param("2 1 -v");
    model = pars.parse();
    CPPUNIT_ASSERT( pars.version() ); 
  }

  void testParseSeeds() {
    Model model;
    Param pars = Param("4 7 -seed 123 -t 5");
    model = pars.parse();
    CPPUNIT_ASSERT_EQUAL( (size_t)123, pars.random_seed() );

    pars = Param("4 7 -seed 1 2 3 -t 5");
    model = pars.parse();
    size_t seed = pars.random_seed();
    pars = Param("4 7 -seed 1 2 3 -t 5");
    model = pars.parse();
    CPPUNIT_ASSERT_EQUAL(seed, pars.random_seed());

    CPPUNIT_ASSERT_THROW(Param("20 10 -seed 1 2 3 4").parse(), std::invalid_argument); 
    CPPUNIT_ASSERT_THROW(Param("20 10 -seed -t 2").parse(), std::invalid_argument); 
  }

  void testParseCommonOptions() {
    Model model = Param("4 7 -t 40.04 -r 1.24 1001 -l 1000").parse();
    CPPUNIT_ASSERT_EQUAL( (size_t)4, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)7, model.loci_number() );
    CPPUNIT_ASSERT( areSame(40.04/(4*model.default_pop_size*1001), model.mutation_rate()) );
    CPPUNIT_ASSERT( areSame(1.24/(4*model.default_pop_size*1000), model.recombination_rate()) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1001, model.loci_length() );
    CPPUNIT_ASSERT( model.has_approximation() );
    CPPUNIT_ASSERT( model.has_window_seq() );
    CPPUNIT_ASSERT_EQUAL( 1000.0, model.window_length_seq() );

    CPPUNIT_ASSERT_THROW( Param("15 10 -I 3 7 8 5").parse(), std::invalid_argument ); 
    CPPUNIT_ASSERT_THROW( Param("15 10 -tv").parse(), std::invalid_argument ); 
    CPPUNIT_ASSERT_THROW( Param("15 10 -t 17 1").parse(), std::invalid_argument ); 

    model = Param("20 10 -t 3.74 -I 3 7 8 5 -T -M 5.0").parse();
    CPPUNIT_ASSERT_EQUAL( (size_t)3, model.population_number() );
    CPPUNIT_ASSERT_EQUAL( (size_t)20, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(4), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(10), (size_t)1 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(17), (size_t)2 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(4), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(17), (double)0.0 );

    model = Param("20 10 -t 3.74 -I 3 7 8 5 5.0 -T").parse();
    CPPUNIT_ASSERT_EQUAL( (size_t)3, model.population_number() );
    CPPUNIT_ASSERT_EQUAL( (size_t)20, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(4), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(10), (size_t)1 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(17), (size_t)2 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(4), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(17), (double)0.0 );
    //std::cout << model << std::endl;
    model.finalize();
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 2)) );

    model = Param("2 1 -t 3.74 -I 2 1 1 5.0").parse();
    CPPUNIT_ASSERT_EQUAL( (size_t)2, model.population_number() );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(0), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(1), (size_t)1 );

    model = Param("23 10 -t 3.74 -I 3 7 8 5 -eI 12.3 2 0 1 -M 5.0").parse();
    CPPUNIT_ASSERT_EQUAL( model.sample_population(20), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(22), (size_t)2 );
    CPPUNIT_ASSERT_EQUAL( 12.3 * 4 * model.default_pop_size, model.sample_time(20) );
    CPPUNIT_ASSERT_EQUAL( 12.3 * 4 * model.default_pop_size, model.sample_time(22) );

    // -N & -eN 
    model = Param("20 10 -t 3.74 -I 3 7 8 5 -N 0.3 -eN 8.2 0.75 -G 1.5 -M 5.0").parse();
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
    model = Param("20 10 -t 3.74 -I 3 7 8 5 -G 1.5 -n 2 0.3 -eN 1.1 0.75 -en 2 3 0.1 -eG 1.5 2 -M 5.0").parse();
    model.finalize();
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
    Model model = Param("20 10 -I 2 10 10 -ma x 5 7 x").parse();
    model.resetTime();
    CPPUNIT_ASSERT( areSame(5.0/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(7.0/(4*model.default_pop_size), model.migration_rate(1, 0)) );

    // -ema
    model = Param("20 10 -I 2 10 10 -ema 1.6 x 5 7 x").parse();
    model.resetTime();
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.6 * 4 * model.default_pop_size, model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame(5.0/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(7.0/(4*model.default_pop_size), model.migration_rate(1, 0)) );

    // -M
    model = Param("30 10 -I 3 10 10 10 0 -M 5").parse();
    model.resetTime();
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 2)) );

    // -eM
    model = Param("30 10 -I 3 10 10 10 0 -eM 1.6 5").parse();
    model.resetTime();
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.6 * 4 * model.default_pop_size,  model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(2.5/(4*model.default_pop_size), model.migration_rate(0, 2)) );

    model = Param("20 1 -I 2 13 7 -m 1 2 1.5").parse();
    //std::cout << model << std::endl;
    CPPUNIT_ASSERT( areSame(0.0, model.migration_rate(1, 0)) );
    CPPUNIT_ASSERT( areSame(1.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );

    model = Param("20 1 -I 2 13 7 -m 1 2 1.5 -m 2 1 0.5").parse();
    //std::cout << model << std::endl;
    CPPUNIT_ASSERT( areSame(1.5/(4*model.default_pop_size), model.migration_rate(0, 1)) );
    CPPUNIT_ASSERT( areSame(0.5/(4*model.default_pop_size), model.migration_rate(1, 0)) );

    model = Param("20 1 -I 2 13 7 -em 1.0 2 1 0.5 -em 2.0 2 1 0.6").parse();
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
    Model model;
    CPPUNIT_ASSERT_NO_THROW(model = Param("20 10 -I 2 10 10 1.5 -es 1.6 2 0.5").parse());
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

    CPPUNIT_ASSERT_THROW(Param("20 10 -I 2 10 10 1.5 -es 1.6 2 0.5 -eM 0.9 5.0").parse(), 
                         std::invalid_argument );
    CPPUNIT_ASSERT_THROW(Param("20 10 -I 2 10 10 1.5 -es 1.6 2 0.5 -eG 0.9 5.0").parse(), 
                         std::invalid_argument );
    CPPUNIT_ASSERT_THROW(Param("20 10 -I 2 10 10 1.5 -es 1.6 2 0.5 -eN 0.9 5.0").parse(), 
                         std::invalid_argument );
    CPPUNIT_ASSERT_THROW(Param("20 10 -I 2 10 10 1.5 -es 1.6 2 0.5 -es 0.9 3 1").parse(), 
                         std::invalid_argument );
  }

  void testParseSplitOptions() {
    // -ej
    Model model;
    CPPUNIT_ASSERT_NO_THROW(model = Param("20 10 -I 2 10 10 -ej 1.6 2 1 -M 1.3").parse());
    model.resetTime();
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.6 * 4 * model.default_pop_size, model.getCurrentTime() );
    CPPUNIT_ASSERT_EQUAL( 1.0, model.single_mig_pop(1, 0) );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.single_mig_pop(0, 1) );

    CPPUNIT_ASSERT_EQUAL( 0.0, model.migration_rate(0, 1) );
    CPPUNIT_ASSERT_EQUAL( 1.3/(4*model.default_pop_size), model.migration_rate(1, 0) );
  }

  void testParseOutputOptions() {
    Model model;
    CPPUNIT_ASSERT_NO_THROW(model = Param("20 10 -T").parse(); );
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 1 );

    CPPUNIT_ASSERT_NO_THROW(model = Param("20 10 -O").parse(); );
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 1 );

    CPPUNIT_ASSERT_NO_THROW(model = Param("20 10 -t 5").parse(); );
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 1 );

    CPPUNIT_ASSERT_NO_THROW(model = Param("20 10 -L").parse());
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 1 );

    CPPUNIT_ASSERT_THROW(Param("20 10 -oSFS").parse(), std::invalid_argument ); 

    CPPUNIT_ASSERT_NO_THROW(model = Param("20 10 -oSFS -t 5").parse());
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 2 );

    CPPUNIT_ASSERT_NO_THROW(model = Param("20 10 -oSFS -t 5 -L").parse());
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 3 );

    CPPUNIT_ASSERT_NO_THROW(model = Param("20 10 -oSFS -t 5 -L -T").parse());
    CPPUNIT_ASSERT( model.summary_statistics_.size() == 4 );
  }

  void testParseGrowthOptions() {
    // -G && -eG
    Model model;
    char *argv[] = { "scrm", "20", "10", "-t", "3.74", "-I", "2", "10", "10", 
      "-G", "2", "-eG", "1", "3", "-M", "5.0" };
    CPPUNIT_ASSERT_NO_THROW( model = Param(16, argv).parse(); );
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
    CPPUNIT_ASSERT_NO_THROW( model = Param(21, argv2).parse(); );
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
    Model model;
    char *argv[] = { "scrm", "20", "10", "-t", "3.74", "-st", "10", "5.1", "-st", "4.5", "3.2"};
    CPPUNIT_ASSERT_NO_THROW( model = Param(11, argv).parse(); );
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
    CPPUNIT_ASSERT_NO_THROW( model = Param(12, argv2).parse(); );
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
    CPPUNIT_ASSERT_NO_THROW( model = Param(12, argv3).parse(); );
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
    CPPUNIT_ASSERT_THROW( Param(11, argv).parse(), std::invalid_argument ); 

    char *argv2[] = { "scrm", "20", "10", "-s", "5"};
    CPPUNIT_ASSERT_THROW( Param(5, argv2).parse(), std::invalid_argument ); 
  }



  void testParseSequenceScaling() {
    Model model;
    char *argv1[] = { "scrm", "4", "7", "-SC", "ms"};
    Param pars = Param(5, argv1);
    CPPUNIT_ASSERT_NO_THROW( model = pars.parse() );
    CPPUNIT_ASSERT( model.getSequenceScaling() == ms );

    char *argv2[] = { "scrm", "4", "7", "-SC", "rel"};
    pars = Param(5, argv2);
    CPPUNIT_ASSERT_NO_THROW( model = pars.parse() );
    CPPUNIT_ASSERT( model.getSequenceScaling() == relative );

    char *argv3[] = { "scrm", "4", "7", "-SC", "abs"};
    pars = Param(5, argv3);
    CPPUNIT_ASSERT_NO_THROW( model = pars.parse() );
    CPPUNIT_ASSERT( model.getSequenceScaling() == absolute );

    char *argv4[] = { "scrm", "4", "7", "-SC", "blub"};
    CPPUNIT_ASSERT_THROW( Param(5, argv4).parse(), std::invalid_argument );
    char *argv5[] = { "scrm", "4", "7", "-SC"};
    CPPUNIT_ASSERT_THROW( Param(4, argv5).parse(), std::invalid_argument );
  }


  void testErrorOnInfintieRecRate() {
    Param pars = Param("4 7 -r 0.5 1");
    CPPUNIT_ASSERT_THROW(model = pars.parse(), std::invalid_argument);

    pars = Param("4 7 -r 0.5 0");
    CPPUNIT_ASSERT_THROW(model = pars.parse(), std::invalid_argument);
  }


  void testScientificNotation() {
    Param pars;
    pars = Param("4 7 -r 1e3 1001");
    CPPUNIT_ASSERT_NO_THROW(model = pars.parse());
    CPPUNIT_ASSERT_EQUAL(1.0/(4*model.default_pop_size), model.recombination_rate());

    pars = Param("4 7 -r 0.5 2e3");
    CPPUNIT_ASSERT_NO_THROW(model = pars.parse());
    CPPUNIT_ASSERT_EQUAL((size_t)2000, model.loci_length());

    pars = Param("4 1e3 -r 0.5 2");
    CPPUNIT_ASSERT_NO_THROW(model = pars.parse());
    CPPUNIT_ASSERT_EQUAL((size_t)1000, model.loci_number());
  }

  void testPrintModel() {
    Param pars;
    pars = Param("4 7 -r 1 1001");
    CPPUNIT_ASSERT_NO_THROW(model = pars.parse());
    CPPUNIT_ASSERT_EQUAL(false, pars.print_model());

    pars = Param("4 7 -r 1 1001 -print-model");
    CPPUNIT_ASSERT_NO_THROW(model = pars.parse());
    CPPUNIT_ASSERT_EQUAL(true, pars.print_model());

    pars = Param("4 7 --print-model");
    CPPUNIT_ASSERT_NO_THROW(model = pars.parse());
    CPPUNIT_ASSERT_EQUAL(true, pars.print_model());
  }

  void testApproximation() {
    Model model = Param("2 2 -r 10 100 -l -1").parse();
    CPPUNIT_ASSERT( !model.has_approximation() );
    CPPUNIT_ASSERT( !model.has_window_rec() );
    CPPUNIT_ASSERT( !model.has_window_seq() );

    model = Param("2 2 -r 10 100 -l 10").parse();
    CPPUNIT_ASSERT(  model.has_approximation() );
    CPPUNIT_ASSERT( !model.has_window_rec() );
    CPPUNIT_ASSERT(  model.has_window_seq() );
    CPPUNIT_ASSERT_EQUAL( 10.0, model.window_length_seq() );

    model = Param("2 2 -r 10 100 -l 10r").parse();
    CPPUNIT_ASSERT(  model.has_approximation() );
    CPPUNIT_ASSERT(  model.has_window_rec() );
    CPPUNIT_ASSERT( !model.has_window_seq() );
    CPPUNIT_ASSERT_EQUAL( (size_t)10, model.window_length_rec() );
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestParam );

/*
 * A sample test case which can be used as a template.
 */
#include <iostream>
#include <cmath>
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <stdexcept>

#include "../../src/model.h"
#include "../../src/forest.h"
#include "../../src/summary_statistics/tmrca.h"

class TestModel : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestModel );

  CPPUNIT_TEST( testSetGetMutationRate );
  CPPUNIT_TEST( testAddChangePositions );
  CPPUNIT_TEST( testSetGetRecombinationRate );
  CPPUNIT_TEST( testGetPopulationSize );
  CPPUNIT_TEST( testAddChangeTime );
  CPPUNIT_TEST( testAddSampleSizes );
  CPPUNIT_TEST( testAddPopulationSizes );
  CPPUNIT_TEST( testAddRelativePopulationSizes );
  CPPUNIT_TEST( testAddGrowthRates );
  CPPUNIT_TEST( testAddMigRates );
  CPPUNIT_TEST( testAddMigRate );
  CPPUNIT_TEST( testDebugConstructor );
  CPPUNIT_TEST( testIncreaseTime );
  CPPUNIT_TEST( testGetNextTime );
  CPPUNIT_TEST( testGetters );
  CPPUNIT_TEST( testHasFixedTimeEvent );
  CPPUNIT_TEST( testCheck );
  CPPUNIT_TEST( testPopSizeAfterGrowth );
  CPPUNIT_TEST( testAddSummaryStatistic );
  CPPUNIT_TEST( testSetLocusLength );
  CPPUNIT_TEST( testAddPopToVectorList );
  CPPUNIT_TEST( testAddPopToMatrixList );

  CPPUNIT_TEST_SUITE_END();

 public:
  void testAddChangeTime() {
    Model model = Model();
    std::vector<double> v1 = std::vector<double>(1, 1),
                        v2 = std::vector<double>(1, 2),
                        v3 = std::vector<double>(1, 3);

    // Check basic adding first time;
    CPPUNIT_ASSERT( model.addChangeTime(0) == 0 );
    CPPUNIT_ASSERT( model.change_times_.size() == 1 );
    CPPUNIT_ASSERT( model.change_times_.at(0) == 0 );
    CPPUNIT_ASSERT( model.pop_sizes_list_.size() == 1 );
    model.pop_sizes_list_[0] = v1;

    // Check adding a time at the end;
    CPPUNIT_ASSERT_EQUAL( (size_t)1, model.addChangeTime(3) );
    CPPUNIT_ASSERT( model.change_times_.size() == 2 );
    CPPUNIT_ASSERT( model.change_times_.at(1) == 3 );
    CPPUNIT_ASSERT( model.pop_sizes_list_[0] == v1 );
    CPPUNIT_ASSERT( model.pop_sizes_list_[1].empty() );
    model.pop_sizes_list_[1] = v3;

    // Check adding a time in the middle;
    CPPUNIT_ASSERT_EQUAL( (size_t)1, model.addChangeTime(2) );
    CPPUNIT_ASSERT( model.change_times_.size() == 3 );
    CPPUNIT_ASSERT( model.change_times_.at(0) == 0 );
    CPPUNIT_ASSERT( model.change_times_.at(1) == 2 );
    CPPUNIT_ASSERT( model.change_times_.at(2) == 3 );
    CPPUNIT_ASSERT( model.pop_sizes_list_[0] == v1 );
    CPPUNIT_ASSERT( model.pop_sizes_list_[1].empty() );
    CPPUNIT_ASSERT( model.pop_sizes_list_[2] == v3 );
    model.pop_sizes_list_[1] = v2;

    CPPUNIT_ASSERT( model.pop_sizes_list_[0] == v1 );
    CPPUNIT_ASSERT( model.pop_sizes_list_[1] == v2 );
    CPPUNIT_ASSERT( model.pop_sizes_list_[2] == v3 );

    // Check that we don't add a time twice
    CPPUNIT_ASSERT_EQUAL( (size_t)0, model.addChangeTime(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, model.addChangeTime(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, model.addChangeTime(2) );
    CPPUNIT_ASSERT_EQUAL( (size_t)3, model.addChangeTime(3) );
    CPPUNIT_ASSERT( model.change_times_.size() == 4 );
    CPPUNIT_ASSERT( model.pop_sizes_list_.size() == 4 );
  }

  void testAddSampleSizes() {
    Model model = Model();
    model.addSampleSizes(0.0, std::vector<size_t>(1,3));
    CPPUNIT_ASSERT_EQUAL( model.sample_size(), (size_t)3 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(2), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(1), (double)0.0 );

    model = Model(0);
    model.set_population_number(3);
    model.addSampleSizes(0.0, std::vector<size_t>(3,1));
    CPPUNIT_ASSERT_EQUAL( model.sample_size(), (size_t)3 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(0), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(1), (size_t)1 );
    CPPUNIT_ASSERT_EQUAL( model.sample_population(2), (size_t)2 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(0), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(1), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( model.sample_time(2), (double)0.0 );

    std::vector<size_t> sample_size;
    sample_size.push_back(0);
    sample_size.push_back(0);
    sample_size.push_back(2);
    model.addSampleSizes(1.0, sample_size);
    CPPUNIT_ASSERT_EQUAL( (size_t)5, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, model.sample_population(3) );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, model.sample_population(4) );
    CPPUNIT_ASSERT_EQUAL( 1.0, model.sample_time(3) );
    CPPUNIT_ASSERT_EQUAL( 1.0, model.sample_time(4) );

    Model model2 = Model();
    std::vector<size_t> sample_sizes;
    sample_sizes.push_back(5);
    sample_sizes.push_back(2);
    model2.addSampleSizes(0.0, sample_sizes);
    CPPUNIT_ASSERT_EQUAL( model2.sample_size(), (size_t)7 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_population(0), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_population(4), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_population(5), (size_t)1 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_population(6), (size_t)1 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_time(0), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_time(4), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_time(6), (double)0.0 );
    sample_sizes.clear();
    sample_sizes.push_back(2);
    sample_sizes.push_back(1);
    model2.addSampleSizes(7.4, sample_sizes);
    CPPUNIT_ASSERT_EQUAL( model2.sample_size(), (size_t)10 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_population(7), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_population(8), (size_t)0 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_population(9), (size_t)1 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_time(0), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_time(6), (double)0.0 );
    CPPUNIT_ASSERT_EQUAL( model2.sample_time(8), (double)7.4 );
  }

  void testAddPopulationSizes() {
    Model model = Model();
    model.set_population_number(2);
    model.addPopulationSizes(1, std::vector<double>(2, 4));
    model.resetTime();
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.0, model.getCurrentTime() );
    CPPUNIT_ASSERT( model.population_size(0) == 4 );

    model.addPopulationSizes(1, std::vector<double>(2, 5), true, false);
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.0 * 4 * model.default_pop_size, model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame(model.population_size(0), 5) );

    model.addPopulationSizes(2, 10, true, false);
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 2.0 * 4 * model.default_pop_size, model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame(model.population_size(0), 10) );

    auto pop_sizes2 = std::vector<double>();
    pop_sizes2.push_back(7);
    pop_sizes2.push_back(6);
    model.addPopulationSizes(0.0, pop_sizes2);
    model.resetTime();
    CPPUNIT_ASSERT( areSame(model.population_size(0), 7) );
    CPPUNIT_ASSERT( areSame(model.population_size(1), 6) );

    CPPUNIT_ASSERT_THROW( model.addPopulationSizes(1, std::vector<double>(1, 5)), std::logic_error );
    CPPUNIT_ASSERT_THROW( model.addPopulationSizes(1, std::vector<double>(3, 5)), std::logic_error );
  }

  void testAddRelativePopulationSizes() {
    Model model = Model();
    model.set_population_number(2);
    model.addPopulationSizes(1, std::vector<double>(2, .5), false, true);
    model.resetTime();
    model.increaseTime();

    CPPUNIT_ASSERT_EQUAL( 1.0, model.getCurrentTime() );
    CPPUNIT_ASSERT_EQUAL( 0.5 * model.default_pop_size, model.population_size(0) );

    model.addPopulationSizes(1, std::vector<double>(2, .4), true, true);
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.0 * 4 * model.default_pop_size, model.getCurrentTime() );
    CPPUNIT_ASSERT_EQUAL( 0.4 * model.default_pop_size, model.population_size(0) );

    model.addPopulationSizes(2, 10, true, true);
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 2.0 * 4 * model.default_pop_size, model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame(10.0 * model.default_pop_size, model.population_size(0)) );

    CPPUNIT_ASSERT_THROW( model.addPopulationSizes(3, 0, true, true), std::invalid_argument);
  }

  void testAddGrowthRates() {
    Model model = Model();
    model.set_population_number(2);

    model.addGrowthRates(1, std::vector<double>(2, 1.5));
    CPPUNIT_ASSERT( model.growth_rates_list_.at(1).at(0) == 1.5 );

    std::vector<double> growth_rates2 = std::vector<double>();
    growth_rates2.push_back(2.5);
    growth_rates2.push_back(3.5);
    model.addGrowthRates(0, growth_rates2);
    CPPUNIT_ASSERT( model.growth_rates_list_.at(0).at(0) == 2.5 );
    CPPUNIT_ASSERT( model.growth_rates_list_.at(0).at(1) == 3.5 );

    CPPUNIT_ASSERT_THROW( model.addGrowthRates(1, std::vector<double>(1, 5)), std::logic_error );
    CPPUNIT_ASSERT_THROW( model.addGrowthRates(1, std::vector<double>(3, 5)), std::logic_error );
  }

  void testGetPopulationSize() {
    Model model = Model(2);
    model.set_population_number(2);
    model.addSymmetricMigration(0.0, 1);
    model.addGrowthRates(1.0, std::vector<double>(2, 1.5));
    model.finalize();
    model.resetTime();

    CPPUNIT_ASSERT( areSame( model.default_pop_size,  model.population_size() ) );
    CPPUNIT_ASSERT( areSame( model.default_pop_size,  model.population_size(0) ) );
    CPPUNIT_ASSERT( areSame( model.default_pop_size,  model.population_size(1) ) );
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.0, model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame( model.default_pop_size,  model.population_size(0) ) );
    CPPUNIT_ASSERT( areSame( model.default_pop_size * std::exp( -1.5 ),  model.population_size(0, 2.0) ) );
  }

  void testAddMigRates() {
    Model model = Model(2);
    model.set_population_number(3);

    std::vector<double> rates;
    for (size_t i = 1; i < 9; ++i) {
      rates.push_back(i);
    }
    CPPUNIT_ASSERT_THROW( model.addMigrationRates(1, rates), std::logic_error );
    rates.push_back(9);
    model.addMigrationRates(1, rates);

    model.resetTime();
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 2.0, model.migration_rate(0,1) );
    CPPUNIT_ASSERT_EQUAL( 3.0, model.migration_rate(0,2) );
    CPPUNIT_ASSERT_EQUAL( 4.0, model.migration_rate(1,0) );
    CPPUNIT_ASSERT_EQUAL( 7.0, model.migration_rate(2,0) );
    CPPUNIT_ASSERT_EQUAL( 6.0, model.migration_rate(1,2) );
  }

  void testAddMigRate() {
    Model model = Model(2);
    model.set_population_number(2);
    model.addMigrationRate(0.0, 0, 1, 0.5);
    model.resetTime();
    CPPUNIT_ASSERT_EQUAL( 0.5, model.migration_rate(0,1) );
    CPPUNIT_ASSERT( std::isnan(model.migration_rate(1,0)) );

    model.addMigrationRate(0.0, 0, 1, 0.7);
    CPPUNIT_ASSERT_EQUAL( 0.7, model.migration_rate(0,1) );
    CPPUNIT_ASSERT( std::isnan(model.migration_rate(1,0)) );

    model.addMigrationRate(0.0, 1, 0, 0.9);
    CPPUNIT_ASSERT_EQUAL( 0.7, model.migration_rate(0,1) );
    CPPUNIT_ASSERT_EQUAL( 0.9, model.migration_rate(1,0) );
  }

  void testDebugConstructor() {
    Model model = Model(7);
    CPPUNIT_ASSERT_EQUAL( (size_t)7, model.sample_size() );
    CPPUNIT_ASSERT( model.growth_rate(0) == 0 );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, model.current_time_idx_ );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, model.current_seq_idx_ );
  }

  void testAddChangePositions() {
    Model model = Model(7);
    CPPUNIT_ASSERT_EQUAL((size_t)1, model.change_position_.size());
    CPPUNIT_ASSERT_EQUAL(0.0, model.change_position_[0]);
    CPPUNIT_ASSERT_EQUAL((size_t)1, model.mutation_rates_.size());
    CPPUNIT_ASSERT_EQUAL((size_t)1, model.recombination_rates_.size());

    // 0.0 is already in => no change
    CPPUNIT_ASSERT_EQUAL((size_t)0, model.addChangePosition(0.0));
    CPPUNIT_ASSERT_EQUAL((size_t)1, model.change_position_.size());
    CPPUNIT_ASSERT_EQUAL(0.0, model.change_position_[0]);
    CPPUNIT_ASSERT_EQUAL((size_t)1, model.mutation_rates_.size());
    CPPUNIT_ASSERT_EQUAL((size_t)1, model.recombination_rates_.size());

    // Really add 1.0 ad the end
    CPPUNIT_ASSERT_EQUAL((size_t)1, model.addChangePosition(1.0));
    CPPUNIT_ASSERT_EQUAL((size_t)2, model.change_position_.size());
    CPPUNIT_ASSERT_EQUAL(0.0, model.change_position_[0]);
    CPPUNIT_ASSERT_EQUAL(1.0, model.change_position_[1]);
    CPPUNIT_ASSERT_EQUAL((size_t)2, model.mutation_rates_.size());
    CPPUNIT_ASSERT_EQUAL((size_t)2, model.recombination_rates_.size());

    // Really add 1.0 again => nothing changes
    CPPUNIT_ASSERT_EQUAL((size_t)1, model.addChangePosition(1.0));
    CPPUNIT_ASSERT_EQUAL((size_t)2, model.change_position_.size());
    CPPUNIT_ASSERT_EQUAL(0.0, model.change_position_[0]);
    CPPUNIT_ASSERT_EQUAL(1.0, model.change_position_[1]);
    CPPUNIT_ASSERT_EQUAL((size_t)2, model.mutation_rates_.size());
    CPPUNIT_ASSERT_EQUAL((size_t)2, model.recombination_rates_.size());

    // Add 0.5 in the middle
    CPPUNIT_ASSERT_EQUAL((size_t)1, model.addChangePosition(0.5));
    CPPUNIT_ASSERT_EQUAL((size_t)3, model.change_position_.size());
    CPPUNIT_ASSERT_EQUAL(0.0, model.change_position_[0]);
    CPPUNIT_ASSERT_EQUAL(0.5, model.change_position_[1]);
    CPPUNIT_ASSERT_EQUAL(1.0, model.change_position_[2]);
    CPPUNIT_ASSERT_EQUAL((size_t)3, model.mutation_rates_.size());
    CPPUNIT_ASSERT_EQUAL((size_t)3, model.recombination_rates_.size());

    model.setLocusLength(10.0);
    CPPUNIT_ASSERT_THROW(model.addChangePosition(-0.5), std::invalid_argument);
    CPPUNIT_ASSERT_THROW(model.addChangePosition(10.5), std::invalid_argument);
  }

  void testIncreaseTime() {
    Model model = Model(7);
    model.addGrowthRates(1.0, std::vector<double>(1, 1.5));
    model.addGrowthRates(2.0, std::vector<double>(1, 1));

    CPPUNIT_ASSERT_EQUAL( (size_t)0, model.current_time_idx_ );
    CPPUNIT_ASSERT_NO_THROW( model.increaseTime() );

    CPPUNIT_ASSERT_EQUAL( (size_t)1, model.current_time_idx_ );
    CPPUNIT_ASSERT_NO_THROW( model.increaseTime() );

    CPPUNIT_ASSERT_EQUAL( (size_t)2, model.current_time_idx_ );

    CPPUNIT_ASSERT_THROW( model.increaseTime(), std::out_of_range );
  }

  void testGetNextTime() {
    Model model = Model(7);
    CPPUNIT_ASSERT( model.getNextTime() == DBL_MAX );

    model.addGrowthRates(1.0, std::vector<double>(1, 1.5));
    model.addGrowthRates(2.0, std::vector<double>(1, 1));

    CPPUNIT_ASSERT_EQUAL( 1.0, model.getNextTime() );
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 2.0, model.getNextTime() );
    model.increaseTime();
    CPPUNIT_ASSERT( model.getNextTime() == DBL_MAX );
  }

  void testGetters() {
    Model model = Model(7);
    model.addGrowthRates(1.0, std::vector<double>(1, 1.5));
    model.addPopulationSizes(2.0, std::vector<double>(1, 5000));
    model.addGrowthRates(3.0, std::vector<double>(1, 1));
    model.addPopulationSizes(4.0, std::vector<double>(1, 10000));

    CPPUNIT_ASSERT_EQUAL( (size_t)7, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate(0) );

    CPPUNIT_ASSERT_NO_THROW( model.increaseTime() ); // 1.0
    CPPUNIT_ASSERT_EQUAL( 1.5, model.growth_rate(0) );

    CPPUNIT_ASSERT_NO_THROW( model.increaseTime() ); // 2.0
    CPPUNIT_ASSERT_EQUAL( 1.5, model.growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 5000.0, model.population_size(0) );

    CPPUNIT_ASSERT_NO_THROW( model.increaseTime() ); // 3.0
    CPPUNIT_ASSERT_EQUAL( 1.0, model.growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 5000.0, model.population_size(0) );

    CPPUNIT_ASSERT_NO_THROW( model.increaseTime() ); // 4.0
    CPPUNIT_ASSERT_EQUAL( 1.0, model.growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 10000.0, model.population_size(0) );
  }

  void testHasFixedTimeEvent() {
    Model model = Model();
    model.set_population_number(3);

    model.addSingleMigrationEvent(10, 1, 0, .5);
    model.resetTime();
    CPPUNIT_ASSERT( ! model.hasFixedTimeEvent(1) );
    CPPUNIT_ASSERT( ! model.hasFixedTimeEvent(10) );
    model.increaseTime();
    CPPUNIT_ASSERT( model.hasFixedTimeEvent(10) );
    CPPUNIT_ASSERT( ! model.hasFixedTimeEvent(1) );
    CPPUNIT_ASSERT( ! model.hasFixedTimeEvent(20) );
  }

  void testSetGetMutationRate() {
    Model model = Model(5);
    model.setLocusLength(10);
    model.setMutationRate(0.001);
    CPPUNIT_ASSERT_EQUAL( 0.001, model.mutation_rate() );

    model.setMutationRate(10, false, false);
    CPPUNIT_ASSERT_EQUAL( 10.0, model.mutation_rate() );

    model.setMutationRate(10, true, false);
    CPPUNIT_ASSERT_EQUAL( 1.0, model.mutation_rate() );

    model.setMutationRate(10, false, true);
    CPPUNIT_ASSERT_EQUAL( 10.0/(4*model.default_pop_size), model.mutation_rate() );

    model.setMutationRate(10, true, true);
    CPPUNIT_ASSERT_EQUAL(10.0/(4*model.default_pop_size*10),
                         model.mutation_rate() );

    // Test setting multiple rates
    model.setMutationRate(10, false, false, 0.0);
    model.setMutationRate(5, false, false, 2.0);
    CPPUNIT_ASSERT_EQUAL(10.0, model.mutation_rate());
    model.increaseSequencePosition();
    CPPUNIT_ASSERT_EQUAL(5.0, model.mutation_rate());

    model.setMutationRate(7.5, false, false, 1.0);
    model.resetSequencePosition();
    CPPUNIT_ASSERT_EQUAL(10.0, model.mutation_rate());
    model.increaseSequencePosition();
    CPPUNIT_ASSERT_EQUAL(7.5, model.mutation_rate());
    model.increaseSequencePosition();
    CPPUNIT_ASSERT_EQUAL(5.0, model.mutation_rate());

    // Test if the vector are filled during finalization
    model.resetSequencePosition();
    model.setRecombinationRate(2.5, false, false, 0.5);
    model.finalize();
    CPPUNIT_ASSERT_EQUAL(10.0, model.mutation_rate());
    model.increaseSequencePosition();
    CPPUNIT_ASSERT_EQUAL(10.0, model.mutation_rate());
    model.increaseSequencePosition();
    CPPUNIT_ASSERT_EQUAL(7.5, model.mutation_rate());
    model.increaseSequencePosition();
    CPPUNIT_ASSERT_EQUAL(5.0, model.mutation_rate());
  }

  void testSetGetRecombinationRate() {
    Model model = Model();
    CPPUNIT_ASSERT( !model.has_recombination() );

    model.setLocusLength(101);
    model.setRecombinationRate(0.001);
    CPPUNIT_ASSERT_EQUAL( 0.001, model.recombination_rate() );
    CPPUNIT_ASSERT_EQUAL( (size_t)101, model.loci_length() );
    CPPUNIT_ASSERT( model.has_recombination() );

    model.setRecombinationRate(0.001, false, false);
    CPPUNIT_ASSERT_EQUAL( 0.001, model.recombination_rate() );
    CPPUNIT_ASSERT_EQUAL( (size_t)101, model.loci_length() );
    CPPUNIT_ASSERT( model.has_recombination() );

    model.setRecombinationRate(0.001, true, false);
    CPPUNIT_ASSERT_EQUAL( 0.00001, model.recombination_rate() );
    CPPUNIT_ASSERT_EQUAL( (size_t)101, model.loci_length() );
    CPPUNIT_ASSERT( model.has_recombination() );

    model.setRecombinationRate(0.001, false, true);
    CPPUNIT_ASSERT_EQUAL( 0.001/(4*model.default_pop_size), model.recombination_rate() );
    CPPUNIT_ASSERT_EQUAL( (size_t)101, model.loci_length() );

    model.setRecombinationRate(0.001, true, true);
    CPPUNIT_ASSERT( areSame(0.00001/(4*model.default_pop_size), model.recombination_rate()) );
    CPPUNIT_ASSERT_EQUAL( (size_t)101, model.loci_length() );

    // Test setting multiple rates
    model.setRecombinationRate(10, false, false, 0.0);
    model.setRecombinationRate(5, false, false, 2.0);
    CPPUNIT_ASSERT_EQUAL(10.0, model.recombination_rate());
    model.increaseSequencePosition();
    CPPUNIT_ASSERT_EQUAL(5.0, model.recombination_rate());

    model.setRecombinationRate(7.5, false, false, 1.0);
    model.resetSequencePosition();
    CPPUNIT_ASSERT_EQUAL(10.0, model.recombination_rate());
    model.increaseSequencePosition();
    CPPUNIT_ASSERT_EQUAL(7.5, model.recombination_rate());
    model.increaseSequencePosition();
    CPPUNIT_ASSERT_EQUAL(5.0, model.recombination_rate());

    // Check for error on invalid change positions
    CPPUNIT_ASSERT_THROW( model.setRecombinationRate(-0.001, 100), std::invalid_argument);
    model.setLocusLength(10);
    CPPUNIT_ASSERT_THROW(model.setRecombinationRate(7.5, false, false, -1.0), std::invalid_argument);
    CPPUNIT_ASSERT_THROW(model.setRecombinationRate(7.5, false, false, 11.0), std::invalid_argument);
  }

  void testCheck() {
    Model model = Model(1);
    CPPUNIT_ASSERT_THROW( model.check(), std::invalid_argument );

    model = Model(2);
    model.set_population_number(2);
    CPPUNIT_ASSERT_THROW( model.check(), std::invalid_argument );
    model.addMigrationRate(20, 0, 1, 0.0);
    CPPUNIT_ASSERT_THROW( model.finalize(), std::invalid_argument );
    model.addMigrationRate(10, 0, 1, 5.0);
    CPPUNIT_ASSERT_NO_THROW( model.finalize() );

    model = Model(3);
    model.set_population_number(2);
    model.addSingleMigrationEvent(1, 0, 1, 1);
    CPPUNIT_ASSERT_NO_THROW( model.check() );
  }

  void testPopSizeAfterGrowth() {
    // Growth only
    Model model = Model(5);
    model.set_population_number(2);
    model.addSymmetricMigration(0, 1.0);
    model.addPopulationSizes(0, 1000);
    std::vector<double> growth_rates;
    growth_rates.push_back(-0.5);
    growth_rates.push_back(0.5);
    model.addGrowthRates(1.0, growth_rates);
    model.addGrowthRates(2.5, 2);
    model.addSingleMigrationEvent(3.5, 1, 0, 0.5);
    model.finalize();

    model.increaseTime();
    model.increaseTime();
    double n_1 = 1000*std::exp(0.5*1.5),
           n_2 = 1000*std::exp(-0.5*1.5);
    CPPUNIT_ASSERT( areSame( n_1, model.population_size(0) ) );
    CPPUNIT_ASSERT( areSame( n_2, model.population_size(1) ) );
    model.increaseTime();
    CPPUNIT_ASSERT( areSame( n_1*std::exp(-2), model.population_size(0) ) );
    CPPUNIT_ASSERT( areSame( n_2*std::exp(-2), model.population_size(1) ) );
    CPPUNIT_ASSERT( areSame( n_1*std::exp(-2*2), model.population_size(0, 4.5) ) );
    CPPUNIT_ASSERT( areSame( n_2*std::exp(-2*2), model.population_size(1, 4.5) ) );

    // Growth with a pop size change in between
    model = Model(5);
    model.set_population_number(2);
    model.addSymmetricMigration(0, 1.0);
    growth_rates.clear();
    growth_rates.push_back(-0.5);
    growth_rates.push_back(0.5);
    model.addGrowthRates(1.0, growth_rates);
    model.addSymmetricMigration(1.25, 2.0);
    model.addPopulationSize(1.5, 0, 500);
    model.addGrowthRate(2.5, 1, 2.0);
    model.addPopulationSize(3.0, 0, 500);
    model.addGrowthRate(3.5, 0, 2.0);
    model.finalize();
    //std::cout << model << std::endl;


    model.resetTime();
    double size0 = model.default_pop_size;
    double size1 = model.default_pop_size;
    CPPUNIT_ASSERT( areSame( size0, model.population_size(0) ) );
    CPPUNIT_ASSERT( areSame( size1, model.population_size(1) ) );

    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.0, model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame( size0, model.population_size(0) ) );
    CPPUNIT_ASSERT( areSame( size1, model.population_size(1) ) );

    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.25, model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame( size0*std::exp(0.5*0.25), model.population_size(0) ) );
    CPPUNIT_ASSERT( areSame( size0*std::exp(-0.5*0.25), model.population_size(1) ) );

    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 1.5, model.getCurrentTime() );
    CPPUNIT_ASSERT( areSame( 500, model.population_size(0) ) );
    CPPUNIT_ASSERT( areSame( model.default_pop_size*std::exp(-0.5*0.5), model.population_size(1) ) );

    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 2.5, model.getCurrentTime() );
    size0 = 500*std::exp(0.5*1.0);
    size1 = model.default_pop_size*std::exp(-0.5*1.5);
    CPPUNIT_ASSERT( areSame( size0, model.population_size(0) ) );
    CPPUNIT_ASSERT( areSame( size1, model.population_size(1) ) );

    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 3.0, model.getCurrentTime() );
    size0  = 500;
    size1 *= exp(-2*0.5);
    CPPUNIT_ASSERT( areSame( size0, model.population_size(0) ) );
    CPPUNIT_ASSERT( areSame( size1, model.population_size(1) ) );

    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 3.5, model.getCurrentTime() );
    size0 *= exp(0.5*0.5);
    size1 *= exp(-2*0.5);
    CPPUNIT_ASSERT( areSame( size0, model.population_size(0) ) );
    CPPUNIT_ASSERT( areSame( size1, model.population_size(1) ) );

    model = Model(5);
    model.set_population_number(2);
    model.addSymmetricMigration(0, 1.0);
    model.addPopulationSize(0, 1, 500);
    model.finalize();
    //std::cout << model << std::endl;
    model.resetTime();
    CPPUNIT_ASSERT( areSame( model.default_pop_size, model.population_size(0) ) );
    CPPUNIT_ASSERT( areSame( 500.0, model.population_size(1) ) );
  }

  void testAddSummaryStatistic() {
    Model model = Model(5);
    CPPUNIT_ASSERT(model.countSummaryStatistics() == 0);
    model.addSummaryStatistic(std::make_shared<TMRCA>());
    CPPUNIT_ASSERT(model.countSummaryStatistics() == 1);
  }

  void testSetLocusLength() {
    Model model = Model(5);
    model.setLocusLength(10);
    CPPUNIT_ASSERT_EQUAL((size_t)10, model.loci_length() );
    model.setMutationRate(5.0, true, false);
    CPPUNIT_ASSERT_EQUAL(0.5, model.mutation_rate());

    model.setLocusLength(2000);
    CPPUNIT_ASSERT_EQUAL((size_t)2000, model.loci_length() );
    CPPUNIT_ASSERT_EQUAL(0.0025, model.mutation_rate());

    model.setLocusLength(5000);
    CPPUNIT_ASSERT_EQUAL((size_t)5000, model.loci_length() );
    CPPUNIT_ASSERT_EQUAL(0.001, model.mutation_rate());
  }

  void testAddPopToVectorList() {
    Model model = Model(5);
    std::vector<std::vector<double> > vector_list;
    vector_list.push_back(std::vector<double>(2, 1.0));
    vector_list.push_back(std::vector<double>());
    vector_list.push_back(std::vector<double>(2, 0.5));

    model.addPopToVectorList(vector_list);
    CPPUNIT_ASSERT_EQUAL( (size_t)3, vector_list.at(0).size() );
    CPPUNIT_ASSERT_EQUAL( 1.0, vector_list.at(0).at(0) );
    CPPUNIT_ASSERT_EQUAL( 1.0, vector_list.at(0).at(1) );
    CPPUNIT_ASSERT( std::isnan(vector_list.at(0).at(2)));

    CPPUNIT_ASSERT_EQUAL( (size_t)3, vector_list.at(2).size() );
    CPPUNIT_ASSERT_EQUAL( 0.5, vector_list.at(2).at(0) );
    CPPUNIT_ASSERT_EQUAL( 0.5, vector_list.at(2).at(1) );
    CPPUNIT_ASSERT( std::isnan(vector_list.at(2).at(2)));
  }

  void testAddPopToMatrixList() {
    Model model = Model(5);
    std::vector<std::vector<double> > vector_list;
    vector_list.push_back(std::vector<double>(2, 1.0));
    vector_list.push_back(std::vector<double>());
    vector_list.push_back(std::vector<double>(2, 0.5));

    model.set_population_number(3);
    model.addPopToMatrixList(vector_list, 2);

    CPPUNIT_ASSERT_EQUAL( (size_t)6, vector_list.at(0).size() );
    CPPUNIT_ASSERT_EQUAL( 1.0, vector_list.at(0).at(0) );
    CPPUNIT_ASSERT_EQUAL( 1.0, vector_list.at(0).at(2) );
    CPPUNIT_ASSERT( std::isnan(vector_list.at(0).at(1)));
    CPPUNIT_ASSERT( std::isnan(vector_list.at(0).at(3)));
    CPPUNIT_ASSERT( std::isnan(vector_list.at(0).at(4)));
    CPPUNIT_ASSERT( std::isnan(vector_list.at(0).at(5)));

    CPPUNIT_ASSERT_EQUAL( (size_t)6, vector_list.at(2).size() );
    CPPUNIT_ASSERT_EQUAL( 0.5, vector_list.at(2).at(0) );
    CPPUNIT_ASSERT_EQUAL( 0.5, vector_list.at(2).at(2) );
    CPPUNIT_ASSERT( std::isnan(vector_list.at(0).at(1)));
    CPPUNIT_ASSERT( std::isnan(vector_list.at(0).at(3)));
    CPPUNIT_ASSERT( std::isnan(vector_list.at(0).at(4)));
    CPPUNIT_ASSERT( std::isnan(vector_list.at(0).at(5)));
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestModel );

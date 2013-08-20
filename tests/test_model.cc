/*
 * A sample test case which can be used as a template.
 */
#include <iostream>
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <stdexcept>
#include "../src/model.h"

class TestModel : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestModel );

  CPPUNIT_TEST( testDeleteParList );
  CPPUNIT_TEST( testAddChangeTime );
  CPPUNIT_TEST( testAddSampleSizes );
  CPPUNIT_TEST( testAddPopulationSizes );
  CPPUNIT_TEST( testAddRelativePopulationSizes );
  CPPUNIT_TEST( testAddGrowthRates );
  CPPUNIT_TEST( testAddMigRates );
  CPPUNIT_TEST( testDebugConstructor );
  CPPUNIT_TEST( testIncreaseTime );
  CPPUNIT_TEST( testGetNextTime );
  CPPUNIT_TEST( testGetters );

  CPPUNIT_TEST_SUITE_END();

 public:
  void testAddChangeTime() {
    Model model = Model();
    std::vector<size_t> *v1 = new std::vector<size_t>(1, 1), 
                        *v2 = new std::vector<size_t>(1, 2), 
                        *v3 = new std::vector<size_t>(1, 3); 

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
    CPPUNIT_ASSERT( model.pop_sizes_list_[1] == NULL ); 
    model.pop_sizes_list_[1] = v3;

    // Check adding a time in the middle;
    CPPUNIT_ASSERT_EQUAL( (size_t)1, model.addChangeTime(2) );
    CPPUNIT_ASSERT( model.change_times_.size() == 3 );
    CPPUNIT_ASSERT( model.change_times_.at(0) == 0 ); 
    CPPUNIT_ASSERT( model.change_times_.at(1) == 2 ); 
    CPPUNIT_ASSERT( model.change_times_.at(2) == 3 ); 
    CPPUNIT_ASSERT( model.pop_sizes_list_[0] == v1 ); 
    CPPUNIT_ASSERT( model.pop_sizes_list_[1] == NULL ); 
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

  void testDeleteParList() {
    Model model = Model();
    std::vector<size_t>* pop_sizes = new std::vector<size_t>(1, 5);
    model.pop_sizes_list_.push_back(pop_sizes);
    model.deleteParList(model.pop_sizes_list_);
    CPPUNIT_ASSERT( model.pop_sizes_list_.size() == 0 );
    //CPPUNIT_ASSERT_EQUAL( (size_t)0, sample_sizes->size() ); //EVIL??
  }

  void testAddSampleSizes() {
    Model model = Model();
    model.addSampleSizes(0.0, std::vector<size_t>(1,3));
    CPPUNIT_ASSERT_EQUAL( model.sample_size(), (size_t)3 ); 
    CPPUNIT_ASSERT_EQUAL( model.sample_population(2), (size_t)0 ); 
    CPPUNIT_ASSERT_EQUAL( model.sample_time(1), (double)0.0 ); 

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
    model.addPopulationSizes(1, std::vector<size_t>(2, 5));
    CPPUNIT_ASSERT( model.pop_sizes_list_.at(1)->at(0) == 5 );

    std::vector<size_t> pop_sizes2 = std::vector<size_t>();
    pop_sizes2.push_back(7);
    pop_sizes2.push_back(4);
    model.addPopulationSizes(0, pop_sizes2);
    CPPUNIT_ASSERT( model.pop_sizes_list_.at(0)->at(0) == 7 );
    CPPUNIT_ASSERT( model.pop_sizes_list_.at(0)->at(1) == 4 );
    
    CPPUNIT_ASSERT_THROW( model.addPopulationSizes(1, std::vector<size_t>(1, 5)), std::logic_error );
    CPPUNIT_ASSERT_THROW( model.addPopulationSizes(1, std::vector<size_t>(3, 5)), std::logic_error );
  }
  
  void testAddRelativePopulationSizes() {
    Model model = Model();
    model.set_population_number(2);
    model.addRelativePopulationSizes(1, std::vector<double>(2, .5));
    CPPUNIT_ASSERT_EQUAL( model.pop_sizes_list_.at(1)->at(0), (size_t)(0.5 * model.default_pop_size) );
  }

  void testAddGrowthRates() {
    Model model = Model();
    model.set_population_number(2);

    model.addGrowthRates(1, std::vector<double>(2, 1.5));
    CPPUNIT_ASSERT( model.growth_rates_list_.at(1)->at(0) == 1.5 );

    std::vector<double> growth_rates2 = std::vector<double>();
    growth_rates2.push_back(2.5);
    growth_rates2.push_back(3.5);
    model.addGrowthRates(0, growth_rates2);
    CPPUNIT_ASSERT( model.growth_rates_list_.at(0)->at(0) == 2.5 );
    CPPUNIT_ASSERT( model.growth_rates_list_.at(0)->at(1) == 3.5 );

    CPPUNIT_ASSERT_THROW( model.addGrowthRates(1, std::vector<double>(1, 5)), std::logic_error );
    CPPUNIT_ASSERT_THROW( model.addGrowthRates(1, std::vector<double>(3, 5)), std::logic_error );
  }

  void testAddMigRates() {
    Model model = Model();
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
    CPPUNIT_ASSERT_EQUAL( 2.0, model.migration_rate(1,0) );
    CPPUNIT_ASSERT_EQUAL( 3.0, model.migration_rate(2,0) );
    CPPUNIT_ASSERT_EQUAL( 4.0, model.migration_rate(0,1) );
    CPPUNIT_ASSERT_EQUAL( 7.0, model.migration_rate(0,2) );
    CPPUNIT_ASSERT_EQUAL( 6.0, model.migration_rate(2,1) );

    CPPUNIT_ASSERT_EQUAL( 11.0, model.total_migration_rate(0) );
    CPPUNIT_ASSERT_EQUAL( 10.0, model.total_migration_rate(1) );
    CPPUNIT_ASSERT_EQUAL(  9.0, model.total_migration_rate(2) );
  }

  void testDebugConstructor() {
    Model model = Model(7);
    CPPUNIT_ASSERT_EQUAL( (size_t)7, model.sample_size() );
    CPPUNIT_ASSERT( model.growth_rate(0) == 0 );
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
    CPPUNIT_ASSERT( model.getNextTime() == FLT_MAX );

    model.addGrowthRates(1.0, std::vector<double>(1, 1.5));
    model.addGrowthRates(2.0, std::vector<double>(1, 1));

    CPPUNIT_ASSERT_EQUAL( 1.0, model.getNextTime() );
    model.increaseTime();
    CPPUNIT_ASSERT_EQUAL( 2.0, model.getNextTime() );
    model.increaseTime();
    CPPUNIT_ASSERT( model.getNextTime() == FLT_MAX );
  }

  void testGetters() {
    Model model = Model(7);
    model.addGrowthRates(1.0, std::vector<double>(1, 1.5));
    model.addPopulationSizes(2.0, std::vector<size_t>(1, 5000));
    model.addGrowthRates(3.0, std::vector<double>(1, 1));
    model.addPopulationSizes(4.0, std::vector<size_t>(1, 10000));

    CPPUNIT_ASSERT_EQUAL( (size_t)7, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( 0.0, model.growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)10000, model.population_size(0) );

    CPPUNIT_ASSERT_NO_THROW( model.increaseTime() ); // 1.0
    CPPUNIT_ASSERT_EQUAL( 1.5, model.growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)10000, model.population_size(0) );

    CPPUNIT_ASSERT_NO_THROW( model.increaseTime() ); // 2.0
    CPPUNIT_ASSERT_EQUAL( 1.5, model.growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)5000, model.population_size(0) );

    CPPUNIT_ASSERT_NO_THROW( model.increaseTime() ); // 3.0
    CPPUNIT_ASSERT_EQUAL( 1.0, model.growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)5000, model.population_size(0) );

    CPPUNIT_ASSERT_NO_THROW( model.increaseTime() ); // 4.0
    CPPUNIT_ASSERT_EQUAL( 1.0, model.growth_rate(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)10000, model.population_size(0) );
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestModel );

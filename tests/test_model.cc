/*
 * A sample test case which can be used as a template.
 */
#include <iostream>
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/model.h"

class TestModel : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestModel );

  CPPUNIT_TEST( testDeleteParList );
  CPPUNIT_TEST( testAddChangeTime );
  CPPUNIT_TEST( testAddSampleSizes );
  CPPUNIT_TEST( testAddPopulationSizes );
  CPPUNIT_TEST( testAddGrowthRates );
  CPPUNIT_TEST( testDebugConstructor );

  CPPUNIT_TEST_SUITE_END();

 public:
  void testAddChangeTime() {
    Model model = Model();
    std::vector<size_t> *v1 = new std::vector<size_t>(), 
                        *v2 = new std::vector<size_t>(), 
                        *v3 = new std::vector<size_t>(); 

    // Check basic adding first time;
    CPPUNIT_ASSERT( model.addChangeTime(1) == 0 );
    CPPUNIT_ASSERT( model.change_times_.size() == 1 );
    CPPUNIT_ASSERT( model.change_times_.at(0) == 1 ); 
    CPPUNIT_ASSERT( model.sample_sizes_list_.size() == 1 );
    CPPUNIT_ASSERT( model.population_sizes_list_.size() == 1 );
    model.sample_sizes_list_[0] = v1;

    // Check adding a time at the end;
    CPPUNIT_ASSERT_EQUAL( (size_t)1, model.addChangeTime(3) );
    CPPUNIT_ASSERT( model.change_times_.size() == 2 );
    CPPUNIT_ASSERT( model.change_times_.at(1) == 3 ); 
    CPPUNIT_ASSERT( model.sample_sizes_list_[0] == v1 ); 
    CPPUNIT_ASSERT( model.sample_sizes_list_[1] == NULL ); 
    model.sample_sizes_list_[1] = v3;

    // Check adding a time in the middle;
    CPPUNIT_ASSERT_EQUAL( (size_t)1, model.addChangeTime(2) );
    CPPUNIT_ASSERT( model.change_times_.size() == 3 );
    CPPUNIT_ASSERT( model.change_times_.at(0) == 1 ); 
    CPPUNIT_ASSERT( model.change_times_.at(1) == 2 ); 
    CPPUNIT_ASSERT( model.change_times_.at(2) == 3 ); 
    CPPUNIT_ASSERT( model.sample_sizes_list_[0] == v1 ); 
    CPPUNIT_ASSERT( model.sample_sizes_list_[1] == NULL ); 
    CPPUNIT_ASSERT( model.sample_sizes_list_[2] == v3 ); 
    model.sample_sizes_list_[1] = v2;

    // Check adding a time in the beginning;
    CPPUNIT_ASSERT_EQUAL( (size_t)0, model.addChangeTime(0) );
    CPPUNIT_ASSERT( model.change_times_.size() == 4 );
    CPPUNIT_ASSERT( model.change_times_.at(0) == 0 ); 
    CPPUNIT_ASSERT( model.change_times_.at(1) == 1 );
    CPPUNIT_ASSERT( model.sample_sizes_list_[0] == NULL ); 
    CPPUNIT_ASSERT( model.sample_sizes_list_[1] == v1 ); 
    CPPUNIT_ASSERT( model.sample_sizes_list_[2] == v2 ); 
    CPPUNIT_ASSERT( model.sample_sizes_list_[3] == v3 ); 

    // Check that we don't add a time twice 
    CPPUNIT_ASSERT_EQUAL( (size_t)0, model.addChangeTime(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, model.addChangeTime(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, model.addChangeTime(2) );
    CPPUNIT_ASSERT_EQUAL( (size_t)3, model.addChangeTime(3) );
    CPPUNIT_ASSERT( model.change_times_.size() == 4 );
    CPPUNIT_ASSERT( model.sample_sizes_list_.size() == 4 );
    CPPUNIT_ASSERT( model.population_sizes_list_.size() == 4 );
  }

  void testDeleteParList() {
    Model model = Model();
    std::vector<size_t>* sample_sizes = new std::vector<size_t>(1, 5);
    model.sample_sizes_list_.push_back(sample_sizes);
    model.deleteParList(model.sample_sizes_list_);
    CPPUNIT_ASSERT( model.sample_sizes_list_.size() == 0 );
    //CPPUNIT_ASSERT_EQUAL( (size_t)0, sample_sizes->size() ); //EVIL??
  }

  void testAddSampleSizes() {
    Model model = Model();
    std::vector<size_t> sample_sizes = std::vector<size_t>(1, 5);
    model.addSampleSizes(1, sample_sizes);
    CPPUNIT_ASSERT( model.sample_sizes_list_.at(0)->at(0) == 5 );
    std::vector<size_t>* sample_sizes2 = new std::vector<size_t>();
    sample_sizes2->push_back(7);
    sample_sizes2->push_back(5);
    model.addSampleSizes(0, *sample_sizes2);
    delete sample_sizes2;
    CPPUNIT_ASSERT( model.sample_sizes_list_.at(0)->at(0) == 7 );
    CPPUNIT_ASSERT( model.sample_sizes_list_.at(0)->at(1) == 5 );
  }

  void testAddPopulationSizes() {
    Model model = Model();
    std::vector<size_t> pop_sizes = std::vector<size_t>(1, 5);
    model.addPopulationSizes(1, pop_sizes);
    CPPUNIT_ASSERT( model.population_sizes_list_.at(0)->at(0) == 5 );
    std::vector<size_t> pop_sizes2 = std::vector<size_t>();
    pop_sizes2.push_back(7);
    pop_sizes2.push_back(4);
    model.addPopulationSizes(0, pop_sizes2);
    CPPUNIT_ASSERT( model.population_sizes_list_.at(0)->at(0) == 7 );
    CPPUNIT_ASSERT( model.population_sizes_list_.at(0)->at(1) == 4 );
  }

  void testAddGrowthRates() {
    Model model = Model();
    std::vector<double> growth_rates = std::vector<double>(1, 1.5);
    model.addGrowthRates(1, growth_rates);
    CPPUNIT_ASSERT( model.growth_rates_list_.at(0)->at(0) == 1.5 );
    std::vector<double> growth_rates2 = std::vector<double>();
    growth_rates2.push_back(2.5);
    growth_rates2.push_back(3.5);
    model.addGrowthRates(0, growth_rates2);
    CPPUNIT_ASSERT( model.growth_rates_list_.at(0)->at(0) == 2.5 );
    CPPUNIT_ASSERT( model.growth_rates_list_.at(0)->at(1) == 3.5 );
  }

  void testDebugConstructor() {
    Model model = Model(7);
    CPPUNIT_ASSERT_EQUAL( (size_t)7, model.sample_size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)10000, model.population_size() );
  }

  void testGettersAndSetters() {
    //model.set_exact_window_length(100);
    // CPPUNIT_ASSERT( model.exact_window_length() == 100 );
  }

  
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestModel );

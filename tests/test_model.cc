/*
 * A sample test case which can be used as a template.
 */
#include <iostream>
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/model.h"

class TestModel : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestModel );

  CPPUNIT_TEST( testAddTimeFrame );
  CPPUNIT_TEST( testSetTime );
  CPPUNIT_TEST( testDebugConstructor );

  CPPUNIT_TEST_SUITE_END();

 public:
  void testAddTimeFrame() {
    Model model = Model();
    std::set<double>::iterator times;
    TimeFramePars tfp = {10, 100, 1, 1};
    model.addTimeFrame(0, tfp);
    CPPUNIT_ASSERT( model.time_frames_.size() == 1 );  
    CPPUNIT_ASSERT_NO_THROW( model.time_frames_[0] );
    CPPUNIT_ASSERT( model.time_frames_[0].sample_size == 10 );
    times = model.change_times(); 
    CPPUNIT_ASSERT( *times == 0 );

    TimeFramePars tfp2 = {15, 100, 1, 1};
    model.addTimeFrame(5.7, tfp2);
    CPPUNIT_ASSERT( model.time_frames_.size() == 2 );  
    CPPUNIT_ASSERT( model.time_frames_[0].sample_size == 10 );
    CPPUNIT_ASSERT( model.time_frames_[5.7].sample_size == 15 );
    times = model.change_times(); 
    CPPUNIT_ASSERT( *(++times) ==  5.7 );

    TimeFramePars tfp3 = {29, 100, 1, 1};
    model.addTimeFrame(1.2, tfp3);
    CPPUNIT_ASSERT( model.time_frames_.size() == 3 );  
    CPPUNIT_ASSERT( model.time_frames_[0].sample_size == 10 );
    CPPUNIT_ASSERT( model.time_frames_[5.7].sample_size == 15 );
    CPPUNIT_ASSERT( model.time_frames_[1.2].sample_size == 29 );
    times = model.change_times();
    CPPUNIT_ASSERT( *times ==  0 );
    CPPUNIT_ASSERT( *(++times) ==  1.2 );
    CPPUNIT_ASSERT( *(++times) ==  5.7 );
  }

  void testSetTime() {
    Model model = Model();
    TimeFramePars tfp = {10, 100, 1, 1};
    model.addTimeFrame(0, tfp);
    model.setTime(0);
    CPPUNIT_ASSERT( model.current_pars_ == &(model.time_frames_[0]) );
    CPPUNIT_ASSERT( model.current_pars_->sample_size == 10 );
    TimeFramePars tfp2 = {15, 100, 1, 1};
    model.addTimeFrame(5.7, tfp2);
    model.setTime(5.7);
    CPPUNIT_ASSERT( model.current_pars_ == &(model.time_frames_[5.7]) );
    CPPUNIT_ASSERT( model.current_pars_->sample_size == 15 );
  }

  void testDebugConstructor() {
    Model model = Model(7);
    CPPUNIT_ASSERT( model.sample_size() == 7 );
  }

  void testGettersAndSetters() {
    //model.set_exact_window_length(100);
    // CPPUNIT_ASSERT( model.exact_window_length() == 100 );
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestModel );

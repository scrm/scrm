#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../../src/forest.h"
#include "../../src/random/constant_generator.h"
#include "../../src/random/mersenne_twister.h"
#include "../../src/event.h"
#include "../../src/summary_statistics/tmrca.h"

class TestForest : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestForest );

  CPPUNIT_TEST( testInitialization );
  CPPUNIT_TEST( testGettersAndSetters );
  CPPUNIT_TEST( testCreateExampleTree );
  CPPUNIT_TEST( testCheckTreeLength );
  CPPUNIT_TEST( testGetFirstNode );
  CPPUNIT_TEST( testSamplePoint );
  CPPUNIT_TEST( testCalcRecombinationRate );
  CPPUNIT_TEST( testCalcRate );
  CPPUNIT_TEST( testNodeIsOld );
  CPPUNIT_TEST( testPrune );
  CPPUNIT_TEST( testSelectFirstTime );
  CPPUNIT_TEST( testSampleEventType );
  CPPUNIT_TEST( testSampleEvent );
  CPPUNIT_TEST( testGetNodeState ); 
  CPPUNIT_TEST( testCut );
  CPPUNIT_TEST( testImplementCoalescence );
  CPPUNIT_TEST( testBuildInitialTree ); 
  CPPUNIT_TEST( testImplementRecombination ); 
  CPPUNIT_TEST( testImplementFixTimeEvent ); 
  CPPUNIT_TEST( testPrintTree );
  CPPUNIT_TEST( testCopyConstructor );
  CPPUNIT_TEST( testCheckForNodeAtHeight );
  CPPUNIT_TEST( testPrintLocusSumStats );
  CPPUNIT_TEST( testSampleNextPosition ); 
  CPPUNIT_TEST( testClear ); 

  CPPUNIT_TEST_SUITE_END();

 private:
  Model *model, *model_2pop;
  Forest *forest, *forest_2pop;
  MersenneTwister *rg;

 public:
  void setUp() {
    rg = new MersenneTwister(1234);
    model = new Model(5);
    model_2pop = new Model(5);
    model_2pop->set_population_number(2);
    forest = new Forest(model, rg);
    forest->createExampleTree();
    forest_2pop = new Forest(model_2pop, rg);
    forest_2pop->createExampleTree();
  }

  void tearDown() {
    delete forest, forest_2pop;
    delete model, model_2pop;
    delete rg;
  }

  void testInitialization() {
    Forest test_forest = Forest(new Model(4), rg);
    CPPUNIT_ASSERT( test_forest.model().sample_size() == 4 );
    CPPUNIT_ASSERT( test_forest.random_generator() == rg );
    delete test_forest.writable_model();
  }

  void testGettersAndSetters() {
    CPPUNIT_ASSERT( forest->model().sample_size() == 4 );
  }

  void testGetFirstNode() {
    CPPUNIT_ASSERT( forest->nodes()->get(0)->height() == 0 );
  }

  void testCreateExampleTree() {
    CPPUNIT_ASSERT_EQUAL((size_t)2, forest_2pop->model().population_number() );
    CPPUNIT_ASSERT_EQUAL((size_t)9, forest->nodes()->size());
    CPPUNIT_ASSERT( forest->local_root() == forest->nodes()->get(8) );
    CPPUNIT_ASSERT( forest->primary_root() == forest->nodes()->get(8) );
    CPPUNIT_ASSERT( forest->getLocalTreeLength() == 24 );
    CPPUNIT_ASSERT( forest->checkTree() );
    
    forest->createExampleTree();
    CPPUNIT_ASSERT_EQUAL((size_t)9, forest->nodes()->size());
    CPPUNIT_ASSERT( forest->local_root() == forest->nodes()->get(8) );
    CPPUNIT_ASSERT( forest->primary_root() == forest->nodes()->get(8) );
    CPPUNIT_ASSERT( forest->getLocalTreeLength() == 24 );
    CPPUNIT_ASSERT( forest->checkTree() );
  }

  void testCheckTreeLength() {
    CPPUNIT_ASSERT( forest->checkTreeLength() );
  }


  void testCalcRate() {
    TimeIntervalIterator tii(forest, forest->nodes()->at(0));
    size_t pop_size = 2*forest->model().population_size(0);
    Node *node1 = new Node(0.1);
    Node *node2 = new Node(0.2);

    forest->set_active_node(0, node1);
    forest->set_active_node(1, node2);
    forest->states_[0] = 0;
    forest->states_[1] = 0;
    forest->recomb_opp_x_within_scrm = 0;
    CPPUNIT_ASSERT_NO_THROW( forest->calcRates(*tii, forest->recomb_opp_x_within_scrm) );
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[0] );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[1] );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[2] );   

    forest->states_[0] = 1;
    CPPUNIT_ASSERT_NO_THROW( forest->calcRates(*tii, forest->recomb_opp_x_within_scrm) );
    CPPUNIT_ASSERT_EQUAL( 4.0/pop_size, forest->rates_[0] );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[1] );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[2] );   

    forest->states_[1] = 1;
    CPPUNIT_ASSERT_NO_THROW( forest->calcRates(*tii, forest->recomb_opp_x_within_scrm) );
    CPPUNIT_ASSERT( areSame(9.0/pop_size, forest->rates_[0]) );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[1] );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[2] );   

    // Coalescence with structure 
    forest_2pop->set_active_node(0, node1);
    forest_2pop->set_active_node(1, node2);
    forest_2pop->states_[0] = 1;
    forest_2pop->states_[1] = 1;
    forest_2pop->recomb_opp_x_within_scrm = 0;

    node1->set_population(1);
    forest_2pop->nodes()->at(1)->set_population(1);
    TimeIntervalIterator tii2(forest_2pop, forest_2pop->nodes()->at(0));
    CPPUNIT_ASSERT_NO_THROW( forest_2pop->calcRates(*tii2, forest_2pop->recomb_opp_x_within_scrm) );
    // Only node2 can coalescence
    CPPUNIT_ASSERT( areSame(4.0/pop_size, forest_2pop->rates_[0]) );
    CPPUNIT_ASSERT_EQUAL( 0.0, forest_2pop->rates_[1] );
    CPPUNIT_ASSERT_EQUAL( 0.0, forest_2pop->rates_[2] );

    std::vector<double> growth(2, 0.0);
    growth.at(1) = 1.0;
    forest_2pop->writable_model()->addGrowthRates(0, growth);
    growth.at(0) = 2.0;
    forest_2pop->writable_model()->addGrowthRates(1, growth);
    TimeIntervalIterator tii3(forest_2pop, forest_2pop->nodes()->at(0));

    CPPUNIT_ASSERT_NO_THROW( forest_2pop->calcRates(*tii3, forest_2pop->recomb_opp_x_within_scrm) );
    CPPUNIT_ASSERT( areSame(3.0/pop_size, forest_2pop->rates_[0]) );   
    CPPUNIT_ASSERT( areSame(1.0/pop_size, forest_2pop->rates_[1]) );   
    CPPUNIT_ASSERT( areSame(0.0/pop_size, forest_2pop->rates_[2]) );   
    CPPUNIT_ASSERT_EQUAL( (size_t)1, forest_2pop->active_nodes_timelines_[0] );   
    CPPUNIT_ASSERT_EQUAL( (size_t)0, forest_2pop->active_nodes_timelines_[1] );   

    forest_2pop->writable_model()->increaseTime();
    CPPUNIT_ASSERT_NO_THROW( forest_2pop->calcRates(*tii3, forest_2pop->recomb_opp_x_within_scrm) );
    //CPPUNIT_ASSERT_EQUAL( 0.0/pop_size, forest_2pop->rates_[0] );   
    //CPPUNIT_ASSERT_EQUAL( 1.0/pop_size, forest_2pop->rates_[1] );   
    //CPPUNIT_ASSERT_EQUAL( 3.0/pop_size, forest_2pop->rates_[2] );   
    CPPUNIT_ASSERT_EQUAL( (size_t)1, forest_2pop->active_nodes_timelines_[0] );   
    CPPUNIT_ASSERT_EQUAL( (size_t)2, forest_2pop->active_nodes_timelines_[1] );   

    node2->set_population(1);
    CPPUNIT_ASSERT_NO_THROW( forest_2pop->calcRates(*tii3, forest_2pop->recomb_opp_x_within_scrm) );
    //CPPUNIT_ASSERT_EQUAL( 0.0/pop_size, forest_2pop->rates_[0] );   
    //CPPUNIT_ASSERT( areSame(3.0/pop_size, forest_2pop->rates_[1]) );   
    //CPPUNIT_ASSERT_EQUAL( 0.0/pop_size, forest_2pop->rates_[2] );   
    CPPUNIT_ASSERT_EQUAL( (size_t)1, forest_2pop->active_nodes_timelines_[0] );   
    CPPUNIT_ASSERT_EQUAL( (size_t)1, forest_2pop->active_nodes_timelines_[1] );   

    delete node1;
    delete node2;
  }

  void testGetNodeState() { 
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(8), 11) == 1 );
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(6), 11) == 2 );
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(7), 11) == 1 );
  }

  void testPrintTree() {
    //CPPUNIT_ASSERT_NO_THROW( forest->printTree() );
  }

  void testSelectFirstTime() {
    double current_time = -1.0;
    size_t time_line = -1;

    forest->selectFirstTime(-1, 1, current_time, time_line);
    CPPUNIT_ASSERT_EQUAL( -1.0, current_time );

    forest->selectFirstTime(5.0, 1, current_time, time_line);
    CPPUNIT_ASSERT_EQUAL( 5.0, current_time );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, time_line );

    forest->selectFirstTime(7.0, 2, current_time, time_line);
    CPPUNIT_ASSERT_EQUAL( 5.0, current_time );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, time_line );

    forest->selectFirstTime(-1, 2, current_time, time_line);
    CPPUNIT_ASSERT_EQUAL( 5.0, current_time );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, time_line );

    forest->selectFirstTime(4.0, 2, current_time, time_line);
    CPPUNIT_ASSERT_EQUAL( 4.0, current_time );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, time_line );
  }

  void testSampleEventType() {
    Event event;
    Forest *forest2 = forest_2pop; 

    forest2->nodes()->at(0)->set_population(1);
    forest2->nodes()->at(1)->set_population(1);
    forest2->writable_model()->addMigrationRates(1, std::vector<double>(4, 5.0), false, true);
    forest2->writable_model()->addGrowthRates(0.2, std::vector<double>(2, 0.000125));
    forest2->writable_model()->addGrowthRates(1, std::vector<double>(2, 0.0));
    forest2->writable_model()->addGrowthRates(2, std::vector<double>(2, 2));
    forest2->writable_model()->finalize();
    forest2->writable_model()->resetTime();

    TimeIntervalIterator tii(forest2, forest2->nodes()->at(0));

    forest2->set_active_node(0, forest2->nodes()->at(0));
    forest2->set_active_node(1, forest2->nodes()->at(2));
    forest2->states_[0] = 1;
    forest2->states_[1] = 0;
    forest2->calcRates(*tii, forest2->recomb_opp_x_within_scrm);
    forest2->active_nodes_timelines_[0] = 0;
    forest2->active_nodes_timelines_[1] = 0;

    forest2->sampleEventType(-1, 0, *tii, event);
    CPPUNIT_ASSERT( event.isNoEvent() );

    forest2->sampleEventType(0.5, 0, *tii, event);
    CPPUNIT_ASSERT_EQUAL( 0.5, event.time() );
    CPPUNIT_ASSERT( event.isCoalescence() );

    // Only coalescence possible
    for (size_t i = 0; i < 1000; ++i) {
      forest2->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() );
      CPPUNIT_ASSERT( event.node() == forest2->active_node(0) );
    };

    // Coalescence of two nodes
    // active_node 0: Pop 1, 2 Contemporaries 
    // active_node 1: Pop 0, 2 Contemporaries
    // => 50% each 
    forest2->states_[1] = 1;
    forest2->calcRates(*tii, forest2->recomb_opp_x_within_scrm);
    size_t count = 0, count2 = 0;
    for (size_t i = 0; i < 10000; ++i) {
      forest2->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() );
      count += (event.node() == forest2->active_node(0));
    };
    CPPUNIT_ASSERT( 4900 < count && count < 5100 ); // ~5000

    // Test with Pw Coalescence
    // active_node 0: Pop 1, 2 Contemporaries 
    // active_node 1: Pop 1, 2 Contemporaries
    // 4/5 Prob of Coalescence, 1/5 of Pw coalescence
    forest2->writable_model()->resetTime();
    forest2->set_active_node(1, forest2->nodes()->at(1));
    forest2->calcRates(*tii, forest2->recomb_opp_x_within_scrm);
    count = 0;
    for (size_t i = 0; i < 10000; ++i) {
      forest2->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() || event.isPwCoalescence() );
      count += event.isCoalescence();
    };
    CPPUNIT_ASSERT( 7900 < count && count < 8100 ); 

    // Test with migration
    // active_node 0: Pop 1, 2 Contemporaries 
    // active_node 1: Pop 1, 2 Contemporaries
    // => Coal: 4/2Ne, Pw Coal: 1/2Ne
    //    Migration: 2 * 5/4Ne
    // => 40% Coal, 10% Pw Coal, 50% Mig 
    forest2->writable_model()->increaseTime();
    forest2->writable_model()->increaseTime();
    forest2->calcRates(*tii, forest2->recomb_opp_x_within_scrm);

    count = 0;
    for (size_t i = 0; i < 10000; ++i) {
      forest2->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() || event.isPwCoalescence() || event.isMigration() );
      count += event.isMigration();
      count2 += event.isCoalescence();
    };
    CPPUNIT_ASSERT( 4900 < count && count < 5100 );
    CPPUNIT_ASSERT( 3900 < count2 && count2 < 4100 );


    // Coalescence and Recombination 
    // active_node 0: Pop 1, 2 Contemporaries => Coal rate: 2 / 2 * Ne = 1/Ne 
    // active_node 1: Pop 1, Recombination    => Rec rate: 10 Bases * 0.4 / 4Ne = 1/Ne    
    // => 50% Recombination, 50% Coalescence 
    forest2->writable_model()->setLocusLength(101);
    forest2->writable_model()->setRecombinationRate(0.4, false, true);
    forest2->set_current_base(10);
    forest2->active_node(1)->make_nonlocal(forest2->current_rec_);
    forest2->set_next_base(20);
    forest2->current_rec_++;
    forest2->states_[0] = 1;
    forest2->states_[1] = 2;
    forest2->writable_model()->resetTime(); // set migration to 0
    forest2->calcRates(*tii, forest2->recomb_opp_x_within_scrm);

    count = 0;
    for (size_t i = 0; i < 100000; ++i) {
      forest2->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() || event.isRecombination() );
      count += event.isRecombination();
    };
    CPPUNIT_ASSERT( 49500 < count && count < 50500 );

    // Other way round
    Node* tmp = forest2->active_node(1);
    forest2->set_active_node(1, forest2->active_node(0));
    forest2->set_active_node(0, tmp);
    forest2->states_[0] = 2;
    forest2->states_[1] = 1;
    forest2->calcRates(*tii, forest2->recomb_opp_x_within_scrm);

    count = 0;
    for (size_t i = 0; i < 100000; ++i) {
      forest2->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() || event.isRecombination() );
      count += event.isRecombination();
    };
    CPPUNIT_ASSERT( 49500 < count && count < 50500 ); // True with >95%


    // Double recombination
    // => 1/3 active_node 0
    //    2/3 active_node 1
    forest2->states_[0] = 2;
    forest2->states_[1] = 2;

    forest2->set_current_base(10.0);
    forest2->set_next_base(15.0);
    forest2->set_next_base(20.0);
    
    forest2->current_rec_ += 2;
    forest2->active_node(0)->make_nonlocal(forest2->current_rec_ - 1);
    forest2->active_node(1)->make_nonlocal(forest2->current_rec_ - 2);
    forest2->calcRates(*tii, forest2->recomb_opp_x_within_scrm);

    count = 0;
    for (size_t i = 0; i < 15000; ++i) {
      forest2->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isRecombination() );
      count += (event.node() == forest2->active_node(0));
    };
    CPPUNIT_ASSERT( 4880 < count && count < 5120 ); // True with >96%


    // Recombination with Rate 0
    // active_node 1: Up to date = rec rate = 0;
    // => always coalescence
    forest2->states_[0] = 1;
    forest2->active_node(0)->make_local();
    forest2->active_node(1)->set_last_update(forest2->current_rec_);
    forest2->calcRates(*tii, forest2->recomb_opp_x_within_scrm);
    for (size_t i = 0; i < 1000; ++i) {
      forest2->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() );
    };

    // Combinations With Growth
    forest2->writable_model()->resetTime(); 
    forest2->writable_model()->increaseTime();
    forest2->writable_model()->increaseTime();
    forest2->writable_model()->increaseTime();
    forest2->states_[0] = 1;
    forest2->states_[1] = 1;
    forest2->active_node(0)->set_population(0);
    forest2->active_node(1)->set_population(1);
    forest2->calcRates(*tii, forest2->recomb_opp_x_within_scrm);

    for (size_t i = 0; i < 100; ++i) {
      CPPUNIT_ASSERT_NO_THROW( forest2->sampleEventType(0.5, 0, *tii, event) );
      CPPUNIT_ASSERT( event.isMigration() );
      CPPUNIT_ASSERT_NO_THROW( forest2->sampleEventType(0.5, 1, *tii, event) );
      CPPUNIT_ASSERT( event.isCoalescence() );
      CPPUNIT_ASSERT( event.node() == forest2->active_node(0) );
      CPPUNIT_ASSERT_NO_THROW( forest2->sampleEventType(0.5, 2, *tii, event) );
      CPPUNIT_ASSERT( event.isCoalescence() );
      CPPUNIT_ASSERT( event.node() == forest2->active_node(1) );
    };

    forest2->active_node(0)->set_population(0);
    forest2->active_node(1)->set_population(0);
    forest2->calcRates(*tii, forest2->recomb_opp_x_within_scrm);
    for (size_t i = 0; i < 100; ++i) {
      CPPUNIT_ASSERT_NO_THROW( forest2->sampleEventType(0.5, 0, *tii, event) );
      CPPUNIT_ASSERT( event.isMigration() );
      CPPUNIT_ASSERT_NO_THROW( forest2->sampleEventType(0.5, 1, *tii, event) );
      CPPUNIT_ASSERT( event.isCoalescence() || event.isPwCoalescence() );
    };

    delete forest2->writable_model();
    delete forest2;
  }

  void testCalcRecombinationRate() {
    forest->writable_model()->setRecombinationRate(1, false, false, 0);
    forest->writable_model()->setRecombinationRate(2, false, false, 10);
    forest->writable_model()->setRecombinationRate(3, false, false, 15);
    forest->writable_model()->setRecombinationRate(4, false, false, 20);
    forest->writable_model()->resetSequencePosition();

    Node* node = forest->nodes()->at(6);

    forest->set_next_base(7.0);
    forest->current_rec_++;
    CPPUNIT_ASSERT_EQUAL( 2.0, forest->calcRecombinationRate(node) );

    forest->writable_model()->increaseSequencePosition();
    forest->set_current_base(12);
    CPPUNIT_ASSERT_EQUAL( 5.0+2*2, forest->calcRecombinationRate(node) );

    forest->writable_model()->increaseSequencePosition();
    forest->set_current_base(17);
    CPPUNIT_ASSERT_EQUAL( 5.0+10.0+6.0, forest->calcRecombinationRate(node) );

    forest->writable_model()->increaseSequencePosition();
    forest->set_current_base(22);
    CPPUNIT_ASSERT_EQUAL( 5.0+10.0+15.0+8.0, forest->calcRecombinationRate(node) );
  }

  void testSampleEvent() {
    Model model = Model(5);
    Forest forest2 = Forest(&model, rg);
    forest2.createScaledExampleTree();
    forest2.writable_model()->finalize();
    forest2.writable_model()->resetTime();

    TimeIntervalIterator tii(&forest2, forest2.nodes()->at(0));

    forest2.set_active_node(0, forest2.nodes()->at(0));
    forest2.set_active_node(1, forest2.nodes()->at(8));
    forest2.states_[0] = 1;
    forest2.states_[1] = 0;
    double recomb_opp_x_within_scrm = 0;
    forest2.calcRates(*tii, recomb_opp_x_within_scrm);
    forest2.active_nodes_timelines_[0] = 0;
    forest2.active_nodes_timelines_[1] = 0;

    Event event;
    double tmp_event_time = 0.0;
    for (size_t i = 0; i < 1000; ++i) {
      forest2.sampleEvent(*tii, tmp_event_time, event); 
      CPPUNIT_ASSERT( event.isNoEvent() || ( 0 <= event.time() && event.time() < forest2.nodes()->at(4)->height() ) );
      CPPUNIT_ASSERT( event.isNoEvent() || event.isCoalescence() );
    }

    ++tii;
    
    forest2.calcRates(*tii, recomb_opp_x_within_scrm);
    for (size_t i = 0; i < 1000; ++i) {
      forest2.sampleEvent(*tii, tmp_event_time, event); 
      CPPUNIT_ASSERT( event.isNoEvent() || ( forest2.nodes()->at(4)->height() <= event.time() && event.time() < forest2.nodes()->at(5)->height() ) );
      CPPUNIT_ASSERT( event.isNoEvent() || event.isCoalescence() );
    }
  }

  void testNodeIsOld() {
    forest->set_current_base(5.0);
    forest->set_next_base(15);
    forest->current_rec_++;

    forest->writable_model()->disable_approximation();
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(0)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(5)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(6)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(7)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(8)) );

    forest->writable_model()->set_window_length_seq(5);
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(0)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(5)) );
    CPPUNIT_ASSERT(  forest->nodeIsOld(forest->nodes()->at(6)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(7)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(8)) );

    forest->writable_model()->set_window_length_seq(9.5);
    CPPUNIT_ASSERT( forest->nodeIsOld(forest->nodes()->at(6)) );
    forest->writable_model()->set_window_length_seq(10.5);
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(6)) );

    forest->createExampleTree();
    forest->writable_model()->set_window_length_rec(2);
    forest->set_next_base(15);
    forest->current_rec_++;
    forest->set_next_base(20);
    forest->current_rec_++;
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(0)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(5)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(6)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(7)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(8)) );

    forest->set_next_base(25);
    forest->current_rec_++;
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(0)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(5)) );
    CPPUNIT_ASSERT(  forest->nodeIsOld(forest->nodes()->at(6)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(7)) );
    CPPUNIT_ASSERT(! forest->nodeIsOld(forest->nodes()->at(8)) );
  }

  void testPrune() {
    // Old node
    forest->set_current_base(5.0);
    forest->set_next_base(15);
    forest->current_rec_++;

    forest->writable_model()->disable_approximation();
    forest->writable_model()->set_window_length_seq(5);
    CPPUNIT_ASSERT( forest->model().has_approximation() );
    CPPUNIT_ASSERT( forest->model().has_window_seq() );
    CPPUNIT_ASSERT( !forest->model().has_window_rec() );

    CPPUNIT_ASSERT(! forest->pruneNodeIfNeeded(forest->nodes()->at(0)) );
    CPPUNIT_ASSERT(! forest->pruneNodeIfNeeded(forest->nodes()->at(1)) );
    CPPUNIT_ASSERT(! forest->pruneNodeIfNeeded(forest->nodes()->at(2)) );
    CPPUNIT_ASSERT(! forest->pruneNodeIfNeeded(forest->nodes()->at(3)) );
    CPPUNIT_ASSERT(! forest->pruneNodeIfNeeded(forest->nodes()->at(4)) );
    CPPUNIT_ASSERT(! forest->pruneNodeIfNeeded(forest->nodes()->at(5)) );
    CPPUNIT_ASSERT(! forest->pruneNodeIfNeeded(forest->nodes()->at(7)) );
    CPPUNIT_ASSERT(! forest->pruneNodeIfNeeded(forest->nodes()->at(8)) );

    CPPUNIT_ASSERT( forest->pruneNodeIfNeeded(forest->nodes()->at(6)) );
    CPPUNIT_ASSERT( forest->nodes()->size() == 8);
    CPPUNIT_ASSERT( forest->checkTree() );

    // Orphaned node
    CPPUNIT_ASSERT( forest->pruneNodeIfNeeded(forest->nodes()->at(6)) );
    CPPUNIT_ASSERT( forest->nodes()->size() == 7);
    CPPUNIT_ASSERT( forest->checkTree() == 1 );

    // In-Between Nodes should be pruned, iff they are of same age
    Node *parent = new Node(20), 
         *inbetween1 = new Node(19), 
         *inbetween2 = new Node(18), 
         *child = new Node(17);

    forest->set_current_base(13);
    inbetween2->make_nonlocal(forest->current_rec_);
    child->make_nonlocal(forest->current_rec_);

    forest->set_next_base(15);
    forest->current_rec_++;
    parent->make_nonlocal(forest->current_rec_);
    inbetween1->make_nonlocal(forest->current_rec_);
    forest->nodes()->add(parent);
    forest->nodes()->add(inbetween1);
    forest->nodes()->add(inbetween2);
    forest->nodes()->add(child);

    parent->set_first_child(inbetween1);
    inbetween1->set_first_child(inbetween2);
    inbetween2->set_first_child(child);
    child->set_parent(inbetween2);
    inbetween2->set_parent(inbetween1);
    inbetween1->set_parent(parent);

    CPPUNIT_ASSERT( !forest->nodeIsOld(inbetween2) );
    CPPUNIT_ASSERT( forest->pruneNodeIfNeeded(inbetween2) );
    CPPUNIT_ASSERT_EQUAL(inbetween1, child->parent());
    CPPUNIT_ASSERT( inbetween1->parent() == parent );
    CPPUNIT_ASSERT( parent->is_root() );
    CPPUNIT_ASSERT( parent->first_child() == inbetween1 );
    CPPUNIT_ASSERT( inbetween1->first_child() == child );
    CPPUNIT_ASSERT( child->countChildren() == 0 );

    forest->nodes()->at(0)->set_parent(NULL);
    CPPUNIT_ASSERT(! forest->pruneNodeIfNeeded(forest->nodes()->at(0)) );
  }

  void testBuildInitialTree() {
    Model model = Model(5);

    for (size_t i = 0; i < 100; ++i) {
      Forest frst = Forest(&model, rg);
      frst.buildInitialTree();
      CPPUNIT_ASSERT_EQUAL(true, frst.checkTree());
      CPPUNIT_ASSERT( frst.current_base() == 0.0 );
      CPPUNIT_ASSERT( frst.next_base() > 0.0 );
    }

    model = Model(2);
    model.set_population_number(2);
    model.addSampleSizes(0.0, std::vector<size_t>(2, 1)); 
    model.addSymmetricMigration(0.0, 5.0, true, true);
    model.finalize();
    Forest frst = Forest(&model, rg);
    CPPUNIT_ASSERT_NO_THROW( frst.buildInitialTree() );
    CPPUNIT_ASSERT_EQUAL( true, frst.checkTree() );
    size_t i = 0;
    for (size_t j = 0; j < 4; ++j) {
      CPPUNIT_ASSERT( frst.nodes()->at(j)->population() == 0 || 
                      frst.nodes()->at(j)->population() == 1 );
      i += frst.nodes()->at(j)->population();
    }
    CPPUNIT_ASSERT_EQUAL( (size_t)1, i );

    model = Model(0);
    model.set_population_number(2);
    std::vector<size_t> sample_size(2, 1);
    sample_size.at(1) = 0;
    model.addSampleSizes(0.0, sample_size); 
    sample_size.at(0) = 0;
    sample_size.at(1) = 1;
    model.addSampleSizes(1.0, sample_size); 
    model.addSymmetricMigration(0.0, 5.0, true, true);
    model.finalize();
    frst = Forest(&model, rg);
    CPPUNIT_ASSERT_NO_THROW( frst.buildInitialTree() );
    CPPUNIT_ASSERT_EQUAL( true, frst.checkTree() );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, frst.nodes()->at(0)->population() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, frst.nodes()->at(1)->population() );
    CPPUNIT_ASSERT_EQUAL( 0.0, frst.nodes()->at(0)->height() );
    CPPUNIT_ASSERT_EQUAL( 1.0, frst.nodes()->at(1)->height() );
  }

  void testCut() {
    Node* base_node = forest->nodes()->at(4);
    Node* new_root = forest->cut(TreePoint(base_node, 3.5, false));

    CPPUNIT_ASSERT_EQUAL((size_t)11, forest->nodes()->size());
    CPPUNIT_ASSERT( new_root->local() );
    CPPUNIT_ASSERT( new_root->is_root() );
    CPPUNIT_ASSERT_EQUAL((size_t)1, new_root->countChildren() );
    CPPUNIT_ASSERT_EQUAL(3.5, new_root->height() );

    CPPUNIT_ASSERT( base_node->parent() == new_root );
    CPPUNIT_ASSERT( base_node->local() );

    Node* single_branch = forest->local_root()->first_child();
    CPPUNIT_ASSERT( !single_branch->local() );
    CPPUNIT_ASSERT_EQUAL( forest->segment_count(), single_branch->last_update() );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, single_branch->countChildren() );
    CPPUNIT_ASSERT_EQUAL( 3.5, single_branch->height() );
  }

  void testImplementRecombination() {
    Node* new_root = forest->cut(TreePoint(forest->nodes()->at(4), 3.5, false));
    TimeIntervalIterator tii(forest, new_root);
    forest->set_active_node(0, new_root);
    forest->set_active_node(1, forest->local_root());

    forest->states_[0] = 2;
    forest->states_[1] = 0;

    Event event = Event(3.7);
    event.setToCoalescence(new_root, 0);
    forest->tmp_event_ = event;

    forest->implementCoalescence(event, tii);
    ++tii;

    forest->set_current_base(100);
    event.set_time(4.5);
    event.setToRecombination(forest->active_node(0), 0);
    CPPUNIT_ASSERT_NO_THROW( forest->implementRecombination(event, tii) );

    CPPUNIT_ASSERT( forest->nodes()->at(9) == forest->active_node(0) || 
                   forest->nodes()->at(10) == forest->active_node(0) );

    CPPUNIT_ASSERT( forest->active_node(0)->local() );
    CPPUNIT_ASSERT( forest->active_node(0)->is_root() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, forest->active_node(0)->countChildren() );

    Node* single_branch = forest->nodes()->at(9);
    if( forest->nodes()->at(9) == forest->active_node(0) ) 
      single_branch = forest->nodes()->at(10);

    CPPUNIT_ASSERT( !single_branch->local() );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, single_branch->countChildren() );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, single_branch->last_update() );
  }

  void testImplementCoalescence() {
    for (size_t i=0; i<100; ++i) {
      forest->createExampleTree(); 

      Node* new_root = forest->cut(TreePoint(forest->nodes()->at(1), 1.5, false));
      TimeIntervalIterator tii(forest, new_root);
      CPPUNIT_ASSERT_EQUAL( (size_t)3, forest->contemporaries_.size(0) );
      forest->set_active_node(0, new_root);
      forest->set_active_node(1, forest->local_root());

      forest->states_[0] = 1;
      forest->states_[1] = 0;

      forest->tmp_event_ = Event(2.0);
      forest->tmp_event_.setToCoalescence(new_root, 0);

      forest->implementCoalescence(forest->tmp_event_, tii);

      CPPUNIT_ASSERT( forest->active_node(0) == new_root );
      CPPUNIT_ASSERT( new_root->parent_height() == 3 ||
                     new_root->parent_height() == 10 );

      if ( !new_root->local() ) {
        CPPUNIT_ASSERT( new_root->countChildren() == 2 );
        CPPUNIT_ASSERT( !(new_root->first_child()->local() && new_root->second_child()->local()) ); 
        Node* child = new_root->first_child();
        if (!new_root->second_child()->local()) child = new_root->second_child();
        CPPUNIT_ASSERT( child->last_update() == 1 );
        CPPUNIT_ASSERT( child->length_below() == 0 );

        CPPUNIT_ASSERT( new_root->last_update() == 1 );
      }
    }
  }

  void testImplementFixTimeEvent() {
    Model* model2 = new Model(5);
    model2->set_population_number(3);
    model2->addSingleMigrationEvent(0.5, 0, 1, 1.0);
    Forest* forest2 = new Forest(model2, rg);
    forest2->createExampleTree();
    Node* new_root = forest2->cut(TreePoint(forest2->nodes()->at(1), 0.5, false));
    forest2->set_active_node(0, new_root);
    forest2->set_active_node(1, forest2->local_root());
    forest2->states_[0] = 1;
    forest2->states_[1] = 0;
    TimeIntervalIterator tii(forest2, new_root);

    CPPUNIT_ASSERT( new_root->population() == 0 );
    forest2->implementFixedTimeEvent(tii);
    forest2->printTree();
    CPPUNIT_ASSERT( new_root->is_root() );
    CPPUNIT_ASSERT( new_root->population() == 1 );

    // Chained events
    new_root->set_population(0);
    model2->addSingleMigrationEvent(0.5, 1, 2, 1.0);
    forest2->implementFixedTimeEvent(tii);
    forest2->printTree();
    CPPUNIT_ASSERT( new_root->is_root() );
    CPPUNIT_ASSERT( new_root->population() == 2 );

    // Circe detection
    model2->addSingleMigrationEvent(0.5, 2, 0, 1.0);
    CPPUNIT_ASSERT_THROW( forest2->implementFixedTimeEvent(tii), 
                          std::logic_error );

    delete forest2, model2;
  }

  void testSamplePoint() {
    rg->set_seed(1234);
    forest->createScaledExampleTree();

    TreePoint point;
    int n0 = 0, n1 = 0, n2 = 0, n3 = 0, 
        n4 = 0, n5 = 0;
    for (int i = 0; i < 240000; ++i) {
      point = forest->samplePoint();
      CPPUNIT_ASSERT( point.base_node() != NULL );
      CPPUNIT_ASSERT( point.base_node()->local() );
      if (point.base_node() == forest->nodes()->at(0)) ++n0;
      if (point.base_node() == forest->nodes()->at(1)) ++n1;
      if (point.base_node() == forest->nodes()->at(2)) ++n2;
      if (point.base_node() == forest->nodes()->at(3)) ++n3;
      if (point.base_node() == forest->nodes()->at(4)) ++n4;
      if (point.base_node() == forest->nodes()->at(5)) ++n5;
    }

    //std::cout << n0 << " " << n1 << " " << n2 << " "
    //          << n3 << " " << n4 << " " << n5 << std::endl;
    CPPUNIT_ASSERT( 29500 <= n0 && n0 <= 30500 ); // expected 30000 
    CPPUNIT_ASSERT( 29500 <= n1 && n1 <= 30500 ); // expected 30000 
    CPPUNIT_ASSERT(  9800 <= n2 && n2 <= 10200 ); // expected 10000 
    CPPUNIT_ASSERT(  9800 <= n3 && n3 <= 10200 ); // expected 10000 
    CPPUNIT_ASSERT( 89000 <= n4 && n4 <= 91000 ); // expected 90000 
    CPPUNIT_ASSERT( 69000 <= n5 && n5 <= 71000 ); // expected 70000 
  }

  void testCopyConstructor() {
    forest->createScaledExampleTree();
    CPPUNIT_ASSERT( forest->coalescence_finished_ == true );
    Forest forest2 = Forest(*forest);

    CPPUNIT_ASSERT_EQUAL( forest->nodes()->size(), forest2.nodes()->size() );
    CPPUNIT_ASSERT( forest2.model_ == forest->model_ );

    for (auto it = forest2.nodes()->iterator(); it.good(); ++it) {
      CPPUNIT_ASSERT( (*it)->label() <= 4 );
      CPPUNIT_ASSERT( (*it)->parent() != NULL || (*it)->first_child() != NULL );
    }

    CPPUNIT_ASSERT( forest2.checkTree() );
  }

  void testCheckForNodeAtHeight() {
    CPPUNIT_ASSERT( forest->checkForNodeAtHeight( 0.0 ) );
    CPPUNIT_ASSERT( forest->checkForNodeAtHeight( 1.0 ) );
    CPPUNIT_ASSERT( forest->checkForNodeAtHeight( 3.0 ) );
    CPPUNIT_ASSERT( forest->checkForNodeAtHeight( 4.0 ) );
    CPPUNIT_ASSERT( forest->checkForNodeAtHeight( 6.0 ) );
    CPPUNIT_ASSERT( forest->checkForNodeAtHeight( 10.0 ) );

    CPPUNIT_ASSERT( !forest->checkForNodeAtHeight( 0.5 ) );
    CPPUNIT_ASSERT( !forest->checkForNodeAtHeight( 1.5 ) );
    CPPUNIT_ASSERT( !forest->checkForNodeAtHeight( 2.0 ) );
    CPPUNIT_ASSERT( !forest->checkForNodeAtHeight( 9.0 ) );
    CPPUNIT_ASSERT( !forest->checkForNodeAtHeight( 20.0 ) );
  }

  void testPrintLocusSumStats() {
    ostringstream output;
    forest->writable_model()->addSummaryStatistic(std::make_shared<TMRCA>());
    forest->set_current_base(0);
    forest->set_next_base(forest->model().loci_length());
    forest->calcSegmentSumStats();
    forest->printLocusSumStats(output);
    CPPUNIT_ASSERT( output.str() != "" );
  }

  void testSampleNextPosition() {
    forest->createScaledExampleTree();
    forest->writable_model()->setRecombinationRate(0.0);
    forest->writable_model()->setRecombinationRate(1.0, false, false, 3);
    forest->sampleNextBase();
    CPPUNIT_ASSERT_EQUAL(3.0, forest->next_base());
    CPPUNIT_ASSERT_EQUAL(1.0, forest->model().recombination_rate());
  }

  void testClear() {
    forest->createScaledExampleTree();
    forest->writable_model()->setRecombinationRate(0.0);
    forest->writable_model()->setRecombinationRate(1.0, false, false, 3);
    forest->writable_model()->addSymmetricMigration(1.0, 0.5);
    forest->writable_model()->increaseSequencePosition();
    forest->writable_model()->increaseTime();
    CPPUNIT_ASSERT_EQUAL(3.0, forest->model().getCurrentSequencePosition());
    CPPUNIT_ASSERT_EQUAL(1.0, forest->model().getCurrentTime());

    forest->clear();
    CPPUNIT_ASSERT_EQUAL((size_t)0, forest->nodes()->size());
    CPPUNIT_ASSERT_EQUAL((size_t)1, forest->rec_bases_.size());
    CPPUNIT_ASSERT_EQUAL((size_t)0, forest->segment_count());
    CPPUNIT_ASSERT_EQUAL(-1.0, forest->current_base());
    CPPUNIT_ASSERT_EQUAL(0.0, forest->model().getCurrentSequencePosition());
    CPPUNIT_ASSERT_EQUAL(0.0, forest->model().getCurrentTime());
  }
};


//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestForest );

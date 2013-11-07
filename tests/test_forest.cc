/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/forest.h"
#include "../src/random/constant_generator.h"
#include "../src/random/mersenne_twister.h"
#include "../src/event.h"

class TestForest : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestForest );

  CPPUNIT_TEST( testInitialization );
  CPPUNIT_TEST( testGettersAndSetters );
  CPPUNIT_TEST( testCreateExampleTree );
  CPPUNIT_TEST( testCheckTreeLength );
  CPPUNIT_TEST( testGetFirstNode );
  CPPUNIT_TEST( testSamplePoint );
  CPPUNIT_TEST( testCalcRate );
  CPPUNIT_TEST( testIsPrunable );
  CPPUNIT_TEST( testPrune );
  CPPUNIT_TEST( testSelectFirstTime );
  CPPUNIT_TEST( testSampleEventType );
  CPPUNIT_TEST( testSampleEvent );
  CPPUNIT_TEST( testBuildInitialTree ); 
  CPPUNIT_TEST( testCoalescenceWithStructure ); 
  CPPUNIT_TEST( testGetNodeState ); 
  CPPUNIT_TEST( testCut );
  CPPUNIT_TEST( testImplementRecombination ); 
  CPPUNIT_TEST( testPrintTree );
  CPPUNIT_TEST( testCopyConstructor );

  CPPUNIT_TEST_SUITE_END();

 private:
  Forest *forest;
  MersenneTwister *rg;

 public:
  void setUp() {
    rg = new MersenneTwister(1234);
    forest = new Forest(new Model(5), rg);
    forest->createExampleTree();
  }

  void tearDown() {
    delete forest->writable_model();
    delete forest;
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
    CPPUNIT_ASSERT( forest->nodes()->size() == 9 );
    CPPUNIT_ASSERT( forest->local_root() == forest->nodes()->get(8) );
    CPPUNIT_ASSERT( forest->primary_root() == forest->nodes()->get(8) );
    CPPUNIT_ASSERT( forest->local_tree_length() == 24 );
    CPPUNIT_ASSERT( forest->checkTree() == 1 );
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
    CPPUNIT_ASSERT_NO_THROW( forest->calcRates(*tii) );
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[0] );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[1] );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[2] );   

    forest->states_[0] = 1;
    CPPUNIT_ASSERT_NO_THROW( forest->calcRates(*tii) );
    CPPUNIT_ASSERT_EQUAL( 4.0/pop_size, forest->rates_[0] );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[1] );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[2] );   

    forest->states_[1] = 1;
    CPPUNIT_ASSERT_NO_THROW( forest->calcRates(*tii) );
    CPPUNIT_ASSERT( areSame(9.0/pop_size, forest->rates_[0]) );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[1] );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[2] );   

    // Coalescence with structure 
    forest->writable_model()->set_population_number(2);
    node1->set_population(1);
    forest->nodes()->at(1)->set_population(1);
    TimeIntervalIterator tii2(forest, forest->nodes()->at(0));
    CPPUNIT_ASSERT_NO_THROW( forest->calcRates(*tii2) );
    // Only node2 can coalescence
    CPPUNIT_ASSERT( areSame(4.0/pop_size, forest->rates_[0]) );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[1] );   
    CPPUNIT_ASSERT_EQUAL( 0.0, forest->rates_[2] );   

    std::vector<double> growth(2, 0.0);
    growth.at(1) = 1.0;
    forest->writable_model()->addGrowthRates(0, growth);
    growth.at(0) = 2.0;
    forest->writable_model()->addGrowthRates(1, growth);
    TimeIntervalIterator tii3(forest, forest->nodes()->at(0));

    CPPUNIT_ASSERT_NO_THROW( forest->calcRates(*tii3) );
    CPPUNIT_ASSERT_EQUAL( 3.0/pop_size, forest->rates_[0] );   
    CPPUNIT_ASSERT_EQUAL( 1.0/pop_size, forest->rates_[1] );   
    CPPUNIT_ASSERT_EQUAL( 0.0/pop_size, forest->rates_[2] );   
    CPPUNIT_ASSERT_EQUAL( (size_t)1, forest->active_nodes_timelines_[0] );   
    CPPUNIT_ASSERT_EQUAL( (size_t)0, forest->active_nodes_timelines_[1] );   
    
    forest->writable_model()->increaseTime();
    CPPUNIT_ASSERT_NO_THROW( forest->calcRates(*tii3) );
    CPPUNIT_ASSERT_EQUAL( 0.0/pop_size, forest->rates_[0] );   
    CPPUNIT_ASSERT_EQUAL( 1.0/pop_size, forest->rates_[1] );   
    CPPUNIT_ASSERT_EQUAL( 3.0/pop_size, forest->rates_[2] );   
    CPPUNIT_ASSERT_EQUAL( (size_t)1, forest->active_nodes_timelines_[0] );   
    CPPUNIT_ASSERT_EQUAL( (size_t)2, forest->active_nodes_timelines_[1] );   

    node2->set_population(1);
    CPPUNIT_ASSERT_NO_THROW( forest->calcRates(*tii3) );
    CPPUNIT_ASSERT_EQUAL( 0.0/pop_size, forest->rates_[0] );   
    CPPUNIT_ASSERT( areSame(3.0/pop_size, forest->rates_[1]) );   
    CPPUNIT_ASSERT_EQUAL( 0.0/pop_size, forest->rates_[2] );   
    CPPUNIT_ASSERT_EQUAL( (size_t)1, forest->active_nodes_timelines_[0] );   
    CPPUNIT_ASSERT_EQUAL( (size_t)1, forest->active_nodes_timelines_[1] );   

    delete node1;
    delete node2;
  }
  
  void testGetNodeState() { 
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(8), 11) == 1 );
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(6), 11) == 2 );
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(7), 11) == 1 );
  }

  void testPrintTree() {
    CPPUNIT_ASSERT_NO_THROW( forest->printTree() );
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
    Forest *forest2 = new Forest(new Model(5), rg);
    forest2->createExampleTree();
    forest2->writable_model()->set_population_number(2);
    forest2->nodes()->at(0)->set_population(1);
    forest2->nodes()->at(1)->set_population(1);
    forest2->writable_model()->addMigrationRates(1, std::vector<double>(4, 5.0), false, true);
    forest2->writable_model()->addGrowthRates(0.2, std::vector<double>(2, 0.000125));
    forest2->writable_model()->addGrowthRates(1, std::vector<double>(2, 0.0));
    forest2->writable_model()->addGrowthRates(2, std::vector<double>(2, 2));
    forest2->writable_model()->finalize();
    forest2->writable_model()->resetTime();

    TimeIntervalIterator tii(forest2, forest2->nodes()->at(0), false);
    
    forest2->set_active_node(0, forest2->nodes()->at(0));
    forest2->set_active_node(1, forest2->nodes()->at(2));
    forest2->states_[0] = 1;
    forest2->states_[1] = 0;
    forest2->calcRates(*tii);
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
    forest2->calcRates(*tii);
    size_t count = 0, count2 = 0;
    for (size_t i = 0; i < 10000; ++i) {
      forest2->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() );
      count += (event.node() == forest2->active_node(0));
    };
    CPPUNIT_ASSERT( 4950 < count && count < 5050 ); // ~5000

    // Test with Pw Coalescence
    // active_node 0: Pop 1, 2 Contemporaries 
    // active_node 1: Pop 1, 2 Contemporaries
    // 4/5 Prob of Coalescence, 1/5 of Pw coalescence
    forest2->writable_model()->resetTime();
    forest2->set_active_node(1, forest2->nodes()->at(1));
    forest2->calcRates(*tii);
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
    forest2->calcRates(*tii);

    count = 0;
    for (size_t i = 0; i < 10000; ++i) {
      forest2->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() || event.isPwCoalescence() || event.isMigration() );
      count += event.isMigration();
      count2 += event.isCoalescence();
    };
    //std::cout << count << " " << count2 << std::endl;
    CPPUNIT_ASSERT( 4900 < count && count < 5100 );
    CPPUNIT_ASSERT( 3900 < count2 && count2 < 4100 );

    // Test coalescence with recombination 
    // active_node 0: Pop 1, 2 Contemporaries => Coal rate: 2 / 2 * Ne = 1/Ne 
    // active_node 1: Pop 1, Recombination    => Rec rate: 10 Bases * 0.4 / 4Ne = 1/Ne    
    // => 50% Recombination, 50% Coalescence 
    forest2->writable_model()->set_recombination_rate(0.4, 101, false, true);
    forest2->set_current_base(20);
    forest2->states_[0] = 1;
    forest2->states_[1] = 2;
    forest2->active_node(1)->make_nonlocal(10);
    forest2->writable_model()->resetTime(); // set migration to 0
    forest2->calcRates(*tii);
    
    count = 0;
    for (size_t i = 0; i < 10000; ++i) {
      forest2->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() || event.isRecombination() );
      count += event.isRecombination();
    };
    //std::cout << count << std::endl;
    CPPUNIT_ASSERT( 4900 < count && count < 5100 );
    
    // Recombination with Rate 0
    // active_node 1: Up to date = rec rate = 0;
    // => always coalescence
    forest2->active_node(1)->set_last_update(20);
    forest2->calcRates(*tii);
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
    forest2->calcRates(*tii);

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
    forest2->calcRates(*tii);
    for (size_t i = 0; i < 100; ++i) {
      CPPUNIT_ASSERT_NO_THROW( forest2->sampleEventType(0.5, 0, *tii, event) );
      CPPUNIT_ASSERT( event.isMigration() );
      CPPUNIT_ASSERT_NO_THROW( forest2->sampleEventType(0.5, 1, *tii, event) );
      CPPUNIT_ASSERT( event.isCoalescence() || event.isPwCoalescence() );
    };

    delete forest2->writable_model();
    delete forest2;
  }

  void testSampleEvent() {
    Model model = Model(5);
    Forest forest2 = Forest(&model, rg);
    forest2.createScaledExampleTree();
    forest2.writable_model()->finalize();
    forest2.writable_model()->resetTime();

    TimeIntervalIterator tii(&forest2, forest2.nodes()->at(0), false);
    
    forest2.set_active_node(0, forest2.nodes()->at(0));
    forest2.set_active_node(1, forest2.nodes()->at(8));
    forest2.states_[0] = 1;
    forest2.states_[1] = 0;
    forest2.calcRates(*tii);
    forest2.active_nodes_timelines_[0] = 0;
    forest2.active_nodes_timelines_[1] = 0;
    
    Event event;
    double tmp_event_time = 0.0;
    size_t tmp_event_line = -1;
    for (size_t i = 0; i < 1000; ++i) {
      forest2.sampleEvent(*tii, tmp_event_time, tmp_event_line, event); 
      CPPUNIT_ASSERT( event.isNoEvent() || ( 0 <= event.time() && event.time() < forest2.nodes()->at(4)->height() ) );
      CPPUNIT_ASSERT( event.isNoEvent() || event.isCoalescence() );
    }

    ++tii;
    forest2.calcRates(*tii);
    for (size_t i = 0; i < 1000; ++i) {
      forest2.sampleEvent(*tii, tmp_event_time, tmp_event_line, event); 
      CPPUNIT_ASSERT( event.isNoEvent() || ( forest2.nodes()->at(4)->height() <= event.time() && event.time() < forest2.nodes()->at(5)->height() ) );
      CPPUNIT_ASSERT( event.isNoEvent() || event.isCoalescence() );
    }
  }

  void testIsPrunable() {
    forest->writable_model()->set_exact_window_length(5);
    forest->set_current_base(15);
    // Node 6 and 7 are old in this setting, but only node 6 is external and can
    // be pruned.
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(0)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(1)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(2)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(3)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(4)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(5)) );
    CPPUNIT_ASSERT(  forest->isPrunable(forest->getNodes()->get(6)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(7)) );
    CPPUNIT_ASSERT( !forest->isPrunable(forest->getNodes()->get(8)) );

    // Orphaned nodes should be pruned
    Node* orphaned = new Node(12);
    orphaned->make_nonlocal(5);
    forest->nodes()->add(orphaned);
    CPPUNIT_ASSERT( forest->isPrunable(orphaned) );

    // Orphaned active not should not occur, but we wont prune them for safety
    // in case they do...
    // Node* orphaned2 = new Node(14, true);
    // forest->nodes()->add(orphaned2);
    // CPPUNIT_ASSERT( !forest->isPrunable(orphaned2) );

    // In-Between Nodes should be pruned, iff they are of same age
    Node *parent = new Node(20), 
         *inbetween1 = new Node(19), 
         *inbetween2 = new Node(18), 
         *child = new Node(17);

    parent->make_nonlocal(15);
    inbetween1->make_nonlocal(15);
    inbetween2->make_nonlocal(13);
    child->make_nonlocal(13);

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
    
    CPPUNIT_ASSERT( !forest->isPrunable(parent) );
    CPPUNIT_ASSERT( !forest->isPrunable(inbetween1) );
    CPPUNIT_ASSERT(  forest->isPrunable(inbetween2) );
    CPPUNIT_ASSERT( !forest->isPrunable(child) );

    // Local in-between nodes
    parent->make_local();
    inbetween1->make_local();
    inbetween2->make_local();
    child->make_local();
    CPPUNIT_ASSERT( !forest->isPrunable(parent) );
    CPPUNIT_ASSERT(  forest->isPrunable(inbetween1) );
    CPPUNIT_ASSERT(  forest->isPrunable(inbetween2) );
    CPPUNIT_ASSERT( !forest->isPrunable(child) );
  }

  void testPrune() {
    // Old node
    forest->writable_model()->set_exact_window_length(5);
    forest->set_current_base(15);
    forest->prune( forest->nodes()->at(6) );
    CPPUNIT_ASSERT( forest->nodes()->size() == 8);
    CPPUNIT_ASSERT( forest->checkTree() );
    // Orphaned node
    forest->prune( forest->nodes()->at(6) );
    CPPUNIT_ASSERT( forest->nodes()->size() == 7);
    CPPUNIT_ASSERT( forest->checkTree() == 1 );

    // In-Between Nodes should be pruned, iff they are of same age
    Node *parent = new Node(20), 
         *inbetween1 = new Node(19), 
         *inbetween2 = new Node(18), 
         *child = new Node(17);

    parent->make_nonlocal(15);
    inbetween1->make_nonlocal(15);
    inbetween2->make_nonlocal(13);
    child->make_nonlocal(13);
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

    forest->prune(inbetween2);
    CPPUNIT_ASSERT( child->parent() == inbetween1 );
    CPPUNIT_ASSERT( inbetween1->parent() == parent );
    CPPUNIT_ASSERT( parent->is_root() );
    CPPUNIT_ASSERT( parent->first_child() == inbetween1 );
    CPPUNIT_ASSERT( inbetween1->first_child() == child );
    CPPUNIT_ASSERT( child->numberOfChildren() == 0 );
  }

  void testBuildInitialTree() {
    Model* model = new Model();
    model->set_population_number(3);
    std::vector<size_t> sample_sizes;
    sample_sizes.push_back(0);
    sample_sizes.push_back(2);
    sample_sizes.push_back(0);
    model->addSampleSizes(0.0, sample_sizes);
    sample_sizes.at(1) = 0;
    sample_sizes.at(2) = 2;
    model->addSampleSizes(1.0, sample_sizes);

    std::vector<double> mig_rates(9, 1.0);
    model->addMigrationRates(2.0, mig_rates);

    Forest frst = Forest(model, rg);
    frst.buildInitialTree();

    //std::cout << endl;
    //frst.printTree_cout();

    CPPUNIT_ASSERT_EQUAL(true, frst.checkTree());

    CPPUNIT_ASSERT_EQUAL(0.0, frst.nodes()->at(0)->height());
    CPPUNIT_ASSERT_EQUAL(0.0, frst.nodes()->at(1)->height());
    CPPUNIT_ASSERT_EQUAL(1.0, frst.nodes()->at(2)->height());
    CPPUNIT_ASSERT_EQUAL(1.0, frst.nodes()->at(3)->height());

    CPPUNIT_ASSERT_EQUAL((size_t)1, frst.nodes()->at(0)->population());
    CPPUNIT_ASSERT_EQUAL((size_t)1, frst.nodes()->at(1)->population());
    CPPUNIT_ASSERT_EQUAL((size_t)2, frst.nodes()->at(2)->population());
    CPPUNIT_ASSERT_EQUAL((size_t)2, frst.nodes()->at(3)->population());
    
    delete model;
  }

  void testCoalescenceWithStructure() {
    Model* model = new Model();
    model->set_population_number(3);
    std::vector<size_t> sample_sizes;
    sample_sizes.push_back(0);
    sample_sizes.push_back(2);
    sample_sizes.push_back(0);
    model->addSampleSizes(0.0, sample_sizes);
    sample_sizes.at(1) = 0;
    sample_sizes.at(2) = 2;
    model->addSampleSizes(1.0, sample_sizes);

    std::vector<double> mig_rates(9, 1.0);
    model->addMigrationRates(2.0, mig_rates);

    Forest frst = Forest(model, rg);
    frst.buildInitialTree();

    for (NodeIterator ni = frst.nodes()->iterator(); ni.good(); ++ni) {
      CPPUNIT_ASSERT( (*ni)->population() != 0 );
    }

    delete model;
  }

  void testCut() {
    Node* base_node = forest->nodes()->at(4);
    Node* new_root = forest->cut(TreePoint(base_node, 3.5, false));

    CPPUNIT_ASSERT_EQUAL((size_t)11, forest->nodes()->size());
    CPPUNIT_ASSERT( new_root->local() );
    CPPUNIT_ASSERT( new_root->is_root() );
    CPPUNIT_ASSERT_EQUAL(1, new_root->numberOfChildren() );
    CPPUNIT_ASSERT_EQUAL(3.5, new_root->height() );

    CPPUNIT_ASSERT( base_node->parent() == new_root );
    CPPUNIT_ASSERT( base_node->local() );

    Node* single_branch = forest->local_root()->first_child();
    CPPUNIT_ASSERT( !single_branch->local() );
    CPPUNIT_ASSERT_EQUAL( forest->current_base(), single_branch->last_update() );
    CPPUNIT_ASSERT_EQUAL( 0, single_branch->numberOfChildren() );
    CPPUNIT_ASSERT_EQUAL( 3.5, single_branch->height() );
  }

  void testImplementRecombination() {
    Node* new_root = forest->cut(TreePoint(forest->nodes()->at(4), 3.5, false));
    TimeIntervalIterator tii(forest, new_root, false);
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
    CPPUNIT_ASSERT_EQUAL( 1, forest->active_node(0)->numberOfChildren() );

    Node* single_branch = forest->nodes()->at(9);
    if( forest->nodes()->at(9) == forest->active_node(0) ) 
      single_branch = forest->nodes()->at(10);

    //forest->printTree_cout();
    CPPUNIT_ASSERT( !single_branch->local() );
    CPPUNIT_ASSERT_EQUAL( 0, single_branch->numberOfChildren() );
    CPPUNIT_ASSERT_EQUAL( 5.0, single_branch->last_update() );
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

    Forest forest2 = Forest(*forest);

    CPPUNIT_ASSERT( forest2.model_ == forest->model_ );
    for (auto it = forest2.nodes()->iterator(); it.good(); ++it) {
      CPPUNIT_ASSERT( (*it)->label() <= 4 );
      CPPUNIT_ASSERT( (*it)->parent() != NULL || (*it)->first_child() != NULL );
    }

    CPPUNIT_ASSERT( forest2.checkTree() );
    CPPUNIT_ASSERT( forest2.checkLeafsOnLocalTree() );
    CPPUNIT_ASSERT( forest2.checkInvariants() );
  }
};



//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestForest );

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
  //CPPUNIT_TEST( testBuildInitialTree );
  //CPPUNIT_TEST( testCoalescenceWithStructure );
  //CPPUNIT_TEST( testGetNodeState );
  CPPUNIT_TEST( testPrintTree );
  CPPUNIT_TEST_SUITE_END();

 private:
  Forest *forest;
  ConstantGenerator *rg;

 public:
  void setUp() {
    rg = new ConstantGenerator;
    forest = new Forest(new Model(5), rg);
    forest->createExampleTree();
  }

  void tearDown() {
    delete forest;
  }

  void testInitialization() {
    Forest test_forest = Forest(new Model(4), rg);
    CPPUNIT_ASSERT( test_forest.model().sample_size() == 4 );
    CPPUNIT_ASSERT( test_forest.random_generator() == rg );
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


  void testSamplePoint() {
    TreePoint tp = forest->samplePoint();
    CPPUNIT_ASSERT( tp.base_node() != NULL );
    CPPUNIT_ASSERT( tp.relative_height() > 0 );
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
  }
  
  void testGetNodeState() { 
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(0), 5) == 0 );
    CPPUNIT_ASSERT( forest->getNodeState(forest->getNodes()->get(8), 5) == 0 );
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
    MersenneTwister* rg2 = new MersenneTwister(3739);
    forest = new Forest(new Model(5), rg2);
    forest->createExampleTree();
    forest->writable_model()->set_population_number(2);
    forest->nodes()->at(0)->set_population(1);
    forest->nodes()->at(1)->set_population(1);
    forest->writable_model()->addMigrationRates(1, std::vector<double>(4, 0.000125));
    forest->writable_model()->addGrowthRates(0.2, std::vector<double>(2, 0.000125));
    forest->writable_model()->addGrowthRates(1, std::vector<double>(2, 0.0));
    forest->writable_model()->addGrowthRates(2, std::vector<double>(2, 2));
    forest->writable_model()->finalize();
    forest->writable_model()->resetTime();

    TimeIntervalIterator tii(forest, forest->nodes()->at(0));
    
    forest->set_active_node(0, forest->nodes()->at(0));
    forest->set_active_node(1, forest->nodes()->at(2));
    forest->states_[0] = 1;
    forest->states_[1] = 0;
    forest->calcRates(*tii);
    forest->active_nodes_timelines_[0] = 0;
    forest->active_nodes_timelines_[1] = 0;

    forest->sampleEventType(-1, 0, *tii, event);
    CPPUNIT_ASSERT( event.isNoEvent() );

    forest->sampleEventType(0.5, 0, *tii, event);
    CPPUNIT_ASSERT_EQUAL( 0.5, event.time() );
    CPPUNIT_ASSERT( event.isCoalescence() );

    for (size_t i = 0; i < 100; ++i) {
      forest->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() );
      CPPUNIT_ASSERT( event.node() == forest->active_node(0) );
    };

    // Test with Coalescence
    forest->states_[1] = 1;
    forest->calcRates(*tii);
    size_t count = 0, count2 = 0;
    for (size_t i = 0; i < 1000; ++i) {
      forest->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() );
      count += (event.node() == forest->active_node(0));
    };
    CPPUNIT_ASSERT( 450 < count && count < 550 ); // ~500

    // Test with Pw Coalescence
    forest->writable_model()->resetTime();
    forest->set_active_node(1, forest->nodes()->at(1));
    forest->calcRates(*tii);
    count = 0;
    for (size_t i = 0; i < 1000; ++i) {
      forest->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() || event.isPwCoalescence() );
      count += event.isCoalescence();
    };
    CPPUNIT_ASSERT( 700 < count && count < 900 ); // ~800

    // Test with migration
    forest->writable_model()->increaseTime();
    forest->writable_model()->increaseTime();
    forest->calcRates(*tii);
    count = 0;
    for (size_t i = 0; i < 1000; ++i) {
      forest->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() || event.isPwCoalescence() || event.isMigration() );
      count += event.isMigration();
      count2 += event.isCoalescence();
    };
    std::cout << count << std::endl;
    CPPUNIT_ASSERT( 450 < count && count < 550 ); // ~500
    CPPUNIT_ASSERT( 350 < count2 && count2 < 450 ); // ~400

    // Test recombination 
    forest->writable_model()->set_recombination_rate(0.00009, 100);
    forest->set_current_base(20);
    forest->states_[0] = 1;
    forest->states_[1] = 2;
    forest->active_node(1)->make_nonlocal(10);
    forest->writable_model()->resetTime(); // set migration to 0
    forest->calcRates(*tii);
    
    count = 0;
    for (size_t i = 0; i < 1000; ++i) {
      forest->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() || event.isRecombination() );
      count += event.isRecombination();
    };
    CPPUNIT_ASSERT( 875 < count && count < 925 ); // ~500
    
    // Recombination with Rate 0
    forest->active_node(1)->set_last_update(20);
    forest->calcRates(*tii);
    for (size_t i = 0; i < 100; ++i) {
      forest->sampleEventType(0.5, 0, *tii, event);
      CPPUNIT_ASSERT( event.isCoalescence() );
    };

    // Combinations With Growth
    forest->writable_model()->resetTime(); 
    forest->writable_model()->increaseTime();
    forest->writable_model()->increaseTime();
    forest->writable_model()->increaseTime();
    forest->states_[0] = 1;
    forest->states_[1] = 1;
    forest->active_node(0)->set_population(0);
    forest->active_node(1)->set_population(1);
    forest->calcRates(*tii);

    for (size_t i = 0; i < 100; ++i) {
      CPPUNIT_ASSERT_NO_THROW( forest->sampleEventType(0.5, 0, *tii, event) );
      CPPUNIT_ASSERT( event.isMigration() );
      CPPUNIT_ASSERT_NO_THROW( forest->sampleEventType(0.5, 1, *tii, event) );
      CPPUNIT_ASSERT( event.isCoalescence() );
      CPPUNIT_ASSERT( event.node() == forest->active_node(0) );
      CPPUNIT_ASSERT_NO_THROW( forest->sampleEventType(0.5, 2, *tii, event) );
      CPPUNIT_ASSERT( event.isCoalescence() );
      CPPUNIT_ASSERT( event.node() == forest->active_node(1) );
    };

    forest->active_node(0)->set_population(0);
    forest->active_node(1)->set_population(0);
    forest->calcRates(*tii);
    for (size_t i = 0; i < 100; ++i) {
      CPPUNIT_ASSERT_NO_THROW( forest->sampleEventType(0.5, 0, *tii, event) );
      CPPUNIT_ASSERT( event.isMigration() );
      CPPUNIT_ASSERT_NO_THROW( forest->sampleEventType(0.5, 1, *tii, event) );
      CPPUNIT_ASSERT( event.isCoalescence() || event.isPwCoalescence() );
    };
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
    Node* orphaned = new Node(12, false);
    forest->nodes()->add(orphaned);
    CPPUNIT_ASSERT( forest->isPrunable(orphaned) );

    // Orphaned active not should not occur, but we wont prune them for safety
    // in case they do...
    // Node* orphaned2 = new Node(14, true);
    // forest->nodes()->add(orphaned2);
    // CPPUNIT_ASSERT( !forest->isPrunable(orphaned2) );

    // In-Between Nodes should be pruned, iff they are of same age
    Node *parent = new Node(20, false, 15), 
         *inbetween1 = new Node(19, false, 15), 
         *inbetween2 = new Node(18, false, 13), 
         *child = new Node(17, false, 13);
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
    Node *parent = new Node(20, false, 15), 
         *inbetween1 = new Node(19, false, 15), 
         *inbetween2 = new Node(18, false, 13), 
         *child = new Node(17, false, 13);
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

    MersenneTwister* rg2 = new MersenneTwister(123);
    Forest frst = Forest(model, rg2);
    frst.buildInitialTree();

    std::cout << endl;
    frst.printTree_cout();

    CPPUNIT_ASSERT_EQUAL(true, frst.checkTree());

    CPPUNIT_ASSERT_EQUAL(0.0, frst.nodes()->at(0)->height());
    CPPUNIT_ASSERT_EQUAL(0.0, frst.nodes()->at(1)->height());
    CPPUNIT_ASSERT_EQUAL(1.0, frst.nodes()->at(2)->height());
    CPPUNIT_ASSERT_EQUAL(1.0, frst.nodes()->at(3)->height());

    CPPUNIT_ASSERT_EQUAL((size_t)1, frst.nodes()->at(0)->population());
    CPPUNIT_ASSERT_EQUAL((size_t)1, frst.nodes()->at(1)->population());
    CPPUNIT_ASSERT_EQUAL((size_t)2, frst.nodes()->at(2)->population());
    CPPUNIT_ASSERT_EQUAL((size_t)2, frst.nodes()->at(3)->population());
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

    MersenneTwister* rg2 = new MersenneTwister(123);
    Forest frst = Forest(model, rg2);
    frst.buildInitialTree();

    for (NodeIterator ni = frst.nodes()->iterator(); ni.good(); ++ni) {
      CPPUNIT_ASSERT( (*ni)->population() != 0 );
    }
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestForest );

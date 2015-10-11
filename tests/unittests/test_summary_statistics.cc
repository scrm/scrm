/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <valarray>
#include <memory>

#include "../../src/forest.h"
#include "../../src/tree_point.h"
#include "../../src/random/mersenne_twister.h"
#include "../../src/summary_statistics/tmrca.h"
#include "../../src/summary_statistics/seg_sites.h"
#include "../../src/summary_statistics/summary_statistic.h"
#include "../../src/summary_statistics/frequency_spectrum.h"
#include "../../src/summary_statistics/oriented_forest.h"
#include "../../src/summary_statistics/newick_tree.h"

class TestSummaryStatistics : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestSummaryStatistics );

  CPPUNIT_TEST( testTMRCA );
  CPPUNIT_TEST( testSegSitesTraversal );
  CPPUNIT_TEST( testSegSitesGetHaplotypes );
  CPPUNIT_TEST( testSegSitesCalculate );
  CPPUNIT_TEST( testSiteFrequencies );
  CPPUNIT_TEST( testOrientedForestGenerateTreeData );
  CPPUNIT_TEST( testOrientedForest );
  CPPUNIT_TEST( testNewickTree );

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

  void testTMRCA() {
    forest->createScaledExampleTree();
    TMRCA tmrca = TMRCA();
    tmrca.calculate(*forest);
    CPPUNIT_ASSERT_EQUAL( 10.0, tmrca.tmrca().at(0) );
    CPPUNIT_ASSERT_EQUAL( 24.0, tmrca.tree_length().at(0) );

    tmrca.calculate(*forest);
    CPPUNIT_ASSERT_EQUAL( 10.0, tmrca.tmrca().at(1) );
    CPPUNIT_ASSERT_EQUAL( 24.0, tmrca.tree_length().at(1) );

    tmrca.clear();
    CPPUNIT_ASSERT_EQUAL( (size_t)0, tmrca.tmrca().size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, tmrca.tree_length().size() );
  }

  void testSegSitesTraversal() {
    SegSites seg_sites = SegSites();

    std::valarray <bool>haplotype(4);
    seg_sites.traversal(forest->nodes()->at(4), haplotype);
    CPPUNIT_ASSERT_EQUAL( true,  haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( true,  haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[3] );

    haplotype *= 0;
    seg_sites.traversal(forest->nodes()->at(5), haplotype);
    CPPUNIT_ASSERT_EQUAL( false, haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( true,  haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( true,  haplotype[3] );

    haplotype *= 0;
    seg_sites.traversal(forest->nodes()->at(8), haplotype);
    CPPUNIT_ASSERT_EQUAL( true, haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[3] );

    haplotype *= 0;
    seg_sites.traversal(forest->nodes()->at(0), haplotype);
    CPPUNIT_ASSERT_EQUAL( false, haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( true,  haplotype[3] );

    haplotype *= 0;
    seg_sites.traversal(forest->nodes()->at(1), haplotype);
    CPPUNIT_ASSERT_EQUAL( false, haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( true,  haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[3] );
  }

  void testSegSitesGetHaplotypes() {
    SegSites seg_sites = SegSites();

    TreePoint mutation = TreePoint(forest->nodes()->at(4), 1, true);
    std::valarray<bool> haplotype = seg_sites.getHaplotypes(mutation, *forest);
    CPPUNIT_ASSERT_EQUAL( (size_t)4, haplotype.size() );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[3] );

    mutation = TreePoint(forest->nodes()->at(5), 0.5, true);
    haplotype = seg_sites.getHaplotypes(mutation, *forest);
    CPPUNIT_ASSERT_EQUAL( (size_t)4, haplotype.size() );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( false, haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[3] );

    mutation = TreePoint(forest->nodes()->at(8), 0.5, true);
    haplotype = seg_sites.getHaplotypes(mutation, *forest);
    CPPUNIT_ASSERT_EQUAL( (size_t)4, haplotype.size() );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( true, haplotype[3] );
  }

  void testSegSitesCalculate() {
    forest->createScaledExampleTree();
    forest->writable_model()->setMutationRate(0.0001);
    forest->set_current_base(0.0);
    forest->set_next_base(15.0);
    SegSites seg_sites = SegSites();
    CPPUNIT_ASSERT_EQUAL( 0.0, seg_sites.position() );

    // Create an example state
    seg_sites.positions_.push_back(0.5);
    seg_sites.positions_.push_back(0.7);
    seg_sites.positions_.push_back(0.8);

    std::valarray<bool> ht = {1, 0, 0, 0};
    seg_sites.haplotypes_.push_back(ht);
    ht[1] = 1;
    seg_sites.haplotypes_.push_back(ht);
    ht = 0;
    ht[1] = 1;
    seg_sites.haplotypes_.push_back(ht);

    // Check status
    CPPUNIT_ASSERT_EQUAL( (size_t)3, seg_sites.countMutations() );
    CPPUNIT_ASSERT_EQUAL( (size_t)4, seg_sites.getHaplotype(0)->size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)4, seg_sites.getHaplotype(1)->size() );
    CPPUNIT_ASSERT_EQUAL( (size_t)4, seg_sites.getHaplotype(2)->size() );

    // Check clear
    seg_sites.clear();
    CPPUNIT_ASSERT_EQUAL( (size_t)0, seg_sites.countMutations() );

    // Check distribution
    seg_sites.calculate(*forest);
    CPPUNIT_ASSERT( seg_sites.countMutations() > 0 );
    int freqs[4] = { 0 };
    int types[2] = { 0 };
    int sum = 0;
    for (auto it = seg_sites.haplotypes_.begin(); it != seg_sites.haplotypes_.end(); ++it) {
      for (size_t i = 0; i < 4; ++i) {
        if ((*it)[i]) ++freqs[i];
        sum += (*it)[i];
      }
      CPPUNIT_ASSERT( sum == 1 || sum == 2 );
      ++types[sum-1];
      sum = 0;
    }

    for (size_t i = 0; i < 4; ++i) {
      CPPUNIT_ASSERT( 0.4 < (double)freqs[i] / seg_sites.countMutations() );
      CPPUNIT_ASSERT( (double)freqs[i] / seg_sites.countMutations() < 0.43 );
    }

    CPPUNIT_ASSERT( (double)types[0] / types[1] > 0.45 );
    CPPUNIT_ASSERT( (double)types[0] / types[1] < 0.55 );

    CPPUNIT_ASSERT_EQUAL( forest->next_base(), seg_sites.position() );

    // Test if the seg_sites clears itself after each locus automatically
    size_t mutation_count = seg_sites.countMutations();
    rg->set_seed(1234);
    seg_sites.calculate(*forest);
    CPPUNIT_ASSERT_EQUAL( mutation_count, seg_sites.countMutations() );
  }

  void testSiteFrequencies() {
    forest->createScaledExampleTree();
    forest->writable_model()->setMutationRate(0.0001);
    forest->set_current_base(0.0);
    forest->set_next_base(10.0);

    std::shared_ptr<SegSites> seg_sites = std::make_shared<SegSites>();
    seg_sites->positions_.push_back(0.5);
    seg_sites->positions_.push_back(0.7);
    seg_sites->positions_.push_back(0.8);

    std::valarray<bool> ht = {1, 0, 0, 0};
    seg_sites->haplotypes_.push_back(ht); // singleton
    ht[1] = 1;
    seg_sites->haplotypes_.push_back(ht); // doubleton
    ht = 0;
    ht[1] = 1;
    seg_sites->haplotypes_.push_back(ht); // singletion
    seg_sites->set_position(forest->next_base());

    FrequencySpectrum sfs(seg_sites, forest->model());
    std::ostringstream output;

    // Check values for example segment
    sfs.calculate(*forest);
    sfs.printLocusOutput(output);
    CPPUNIT_ASSERT( output.str().compare("SFS: 2 1 0 \n") == 0 );

    // Add another segment
    seg_sites->positions_.push_back(0.9);
    seg_sites->haplotypes_.push_back(ht);
    sfs.calculate(*forest);
    CPPUNIT_ASSERT_EQUAL((size_t)3, sfs.sfs().at(0));
    CPPUNIT_ASSERT_EQUAL((size_t)1, sfs.sfs().at(1));
    CPPUNIT_ASSERT_EQUAL((size_t)0, sfs.sfs().at(2));

    // Check that it is reset at a new locus
    output.str("");
    output.clear();

    sfs.clear();
    forest->createScaledExampleTree();
    forest->set_current_base(0.0);
    forest->set_next_base(0.000001);
    sfs.calculate(*forest);
    sfs.printLocusOutput(output);
    CPPUNIT_ASSERT( output.str().compare("SFS: 0 0 0 \n") == 0 );
  }

  void testOrientedForestGenerateTreeData() {
    OrientedForest of(4);
    size_t pos = 2*forest->sample_size()-2;
    double sf = forest->model().scaling_factor();
    of.generateTreeData(forest->local_root(), pos, 0, sf);
    CPPUNIT_ASSERT( of.heights_.at(0) == 0.0 );
    CPPUNIT_ASSERT( of.heights_.at(1) == 0.0 );
    CPPUNIT_ASSERT( of.heights_.at(2) == 0.0 );
    CPPUNIT_ASSERT( of.heights_.at(3) == 0.0 );
    CPPUNIT_ASSERT( of.parents_.at(4) == 7 );
    CPPUNIT_ASSERT( of.parents_.at(5) == 7 );
    CPPUNIT_ASSERT( of.heights_.at(6) == 10.0 * sf && of.parents_.at(6) == 0 );

    CPPUNIT_ASSERT( of.heights_.at(4) == 1.0 * sf || of.heights_.at(4) == 3.0 * sf );
    if ( of.heights_.at(4) == 1.0 ) {
      CPPUNIT_ASSERT( of.parents_.at(0) == 5 );
      CPPUNIT_ASSERT( of.parents_.at(1) == 5 );
      CPPUNIT_ASSERT( of.parents_.at(2) == 6 );
      CPPUNIT_ASSERT( of.parents_.at(3) == 6 );
      CPPUNIT_ASSERT( of.heights_.at(5) == 3.0 * sf);
    } else {
      CPPUNIT_ASSERT( of.parents_.at(0) == 6 );
      CPPUNIT_ASSERT( of.parents_.at(1) == 6 ); CPPUNIT_ASSERT( of.parents_.at(2) == 5 );
      CPPUNIT_ASSERT( of.parents_.at(3) == 5 );
      CPPUNIT_ASSERT( of.heights_.at(5) == 1.0 * sf );
    }
  }

  void testOrientedForest() {
    forest->createScaledExampleTree();
    forest->set_current_base(0.0);
    forest->set_next_base(10.0);

    OrientedForest of(4);
    std::ostringstream output;
    of.calculate(*forest);
    of.printSegmentOutput(output);
    //std::cout << output.str() << std::endl;
    CPPUNIT_ASSERT( output.str().compare("{\"parents\":[5,5,6,6,7,7,0], \"node_times\":[0,0,0,0,1,3,10]}\n") == 0 ||
                    output.str().compare("{\"parents\":[6,6,5,5,7,7,0], \"node_times\":[0,0,0,0,3,1,10]}\n") == 0 );

    output.str("");
    output.clear();
    forest->writable_model()->setRecombinationRate(0.0001);
    of.clear();
    of.calculate(*forest);
    of.printSegmentOutput(output);
    //std::cout << output.str() << std::endl;
    CPPUNIT_ASSERT( output.str().compare("{\"length\":10, \"parents\":[5,5,6,6,7,7,0], \"node_times\":[0,0,0,0,1,3,10]}\n") == 0 ||
                    output.str().compare("{\"length\":10, \"parents\":[6,6,5,5,7,7,0], \"node_times\":[0,0,0,0,3,1,10]}\n") == 0 );
  }

  void testNewickTree() {
    forest->createScaledExampleTree();
    forest->set_current_base(0.0);
    forest->set_next_base(10.0);

    NewickTree of(4, forest->model().has_recombination());
    std::ostringstream output;
    of.calculate(*forest);
    of.printSegmentOutput(output);
    //std::cout << output.str() << std::endl;
    CPPUNIT_ASSERT( output.str().compare("((1:1,2:1):9,(3:3,4:3):7);\n") == 0 );

    output.str("");
    output.clear();
    forest->writable_model()->setRecombinationRate(0.0001);
    of = NewickTree(4, forest->model().has_recombination());
    of.calculate(*forest);
    of.printSegmentOutput(output);
    //std::cout << output.str() << std::endl;
    CPPUNIT_ASSERT( output.str().compare("[10]((1:1,2:1):9,(3:3,4:3):7);\n") == 0 );

    /** Scientific notation is disabled when using absolute scaling **/
    output.str("");
    output.clear();
    forest->createScaledExampleTree();
    forest->set_current_base(0.0);
    forest->set_next_base(10000000.0);
    of.calculate(*forest);
    of.printSegmentOutput(output);
    //std::cout << output.str() << std::endl;
    CPPUNIT_ASSERT( output.str().compare("[10000000]((1:1,2:1):9,(3:3,4:3):7);\n") == 0 );

    /** Relative sequence scaling **/
    forest->writable_model()->setSequenceScaling(relative);
    forest->writable_model()->setLocusLength(10000);
    forest->createScaledExampleTree();
    forest->set_current_base(0.0);
    forest->set_next_base(100);
    output.str("");
    output.clear();
    of.calculate(*forest);
    of.printSegmentOutput(output);
    //std::cout << output.str() << std::endl;
    CPPUNIT_ASSERT( output.str().compare("[0.01]((1:1,2:1):9,(3:3,4:3):7);\n") == 0 );
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestSummaryStatistics );

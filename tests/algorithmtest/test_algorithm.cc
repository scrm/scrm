/*
 * Long running tests that ensure that scrm produces correct results
 * in various settings.
 */

#pragma GCC diagnostic ignored "-Wwrite-strings"

#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include <valarray>
#include <cmath>

#include "../../src/forest.h"
#include "../../src/param.h"
#include "../../src/random/mersenne_twister.h"


class TestAlgorithm : public CppUnit::TestCase {

	CPPUNIT_TEST_SUITE( TestAlgorithm );

	CPPUNIT_TEST( testInitialTree );
	CPPUNIT_TEST( testARG );
	CPPUNIT_TEST( testPruning );
	CPPUNIT_TEST( testMigration );
	CPPUNIT_TEST( testSizeChange );
	CPPUNIT_TEST( testGrowth );
	CPPUNIT_TEST( testSplit );
	CPPUNIT_TEST( testMerge );

	CPPUNIT_TEST_SUITE_END();

 private:
  MersenneTwister *rg;

  void testTree(Model &model, size_t replicates, 
                double tmrca_mean, double tmrca_sd,
                double tree_length_mean, double tree_length_sd) {
    double tmrca[6] = { 0 }, tree_length[6] = { 0 };
    
    std::cout << "." << std::flush;

    Forest forest = Forest(&model, this->rg);
    for (size_t i = 0; i < replicates; ++i) {
      // Check initial tree
      forest.buildInitialTree();
      tmrca[0] += forest.getTMRCA(true);
      tree_length[0] += forest.getLocalTreeLength(true);

      // Check after recombinations
      if (model.recombination_rate() > 0) {
        for (size_t j = 1; j <= 5; j++) {
          while (forest.next_base() < j*5) {
            forest.sampleNextGenealogy();
          }
          tmrca[j] += forest.getTMRCA(true);
          tree_length[j] += forest.getLocalTreeLength(true);
        }
      }

      // Clear Forest
      forest.clear();
    }

    // Allow an relative error of 2.5%. It would be nice to calculate
    // standard errors, but there's nothing in the std and I'm to lazy to
    // implement very it myself...
    for (int i = 0; i <= 5; ++i) {
      if (i > 0 && tmrca[i] == 0 && tree_length[i] == 0) continue;
      tmrca[i] /= replicates;
      double SE = tmrca_sd / sqrt(replicates);
      if (tmrca[i] < tmrca_mean - 4 * SE || tmrca_mean + 4 * SE < tmrca[i]) {
        std::cout << std::endl 
                  << "TMRCA outside expected range. Observed: " << tmrca[i] 
                  << " Expected: " << tmrca_mean << std::endl;  
        CPPUNIT_ASSERT( false );
      }

      tree_length[i] /= replicates; 
      SE = tree_length_sd / sqrt(replicates);
      if (tree_length[i] < tree_length_mean - 4 * SE || tree_length_mean + 4 * SE < tree_length[i]) {
        std::cout << std::endl 
                  << "Local Tree Length outside expected range. Observed: " << tree_length[i] 
                  << " Expected: " << tree_length_mean << std::endl;  
        CPPUNIT_ASSERT( false );
      }
    }
  }

 public:
  void setUp() {
    rg = new MersenneTwister(78361);
  }

  void tearDown() {
    delete rg;
  }

  void testInitialTree() {
    Model model = Model(5);
    model.setRecombinationRate(0);
    testTree(model, 10000, 0.8, 0.53, 2.08, 1.19); 

    model = Model(10);
    model.setRecombinationRate(0);
    testTree(model, 10000, 0.9, 0.53, 2.83, 1.24); 

    model = Model(20);
    model.setRecombinationRate(0);
    testTree(model, 10000, 0.95, 0.54, 3.55, 1.26); 
  }

  void testARG() {
    Model model = Model(5);
    model.setRecombinationRate(1, false, true);
    testTree(model, 10000, 0.8, 0.53, 2.08, 1.19); 

    model = Model(10);
    model.setRecombinationRate(1, false, true);
    testTree(model, 10000, 0.9, 0.53, 2.83, 1.24); 

    model = Model(20);
    model.setRecombinationRate(1, false, true);
    testTree(model, 10000, 0.95, 0.54, 3.55, 1.26); 
  }

  void testPruning() {
    Model model = Model(10);
    model.setRecombinationRate(1, false, true);
    model.set_window_length_seq(0);
    testTree(model, 10000, 0.9, 0.53, 2.83, 1.24); 

    model.set_window_length_seq(5);
    testTree(model, 5000, 0.9, 0.53, 2.83, 1.24); 

    model.set_window_length_seq(10);
    testTree(model, 5000, 0.9, 0.53, 2.83, 1.24); 

    model.disable_approximation();
    model.set_window_length_rec(0);
    testTree(model, 10000, 0.9, 0.53, 2.83, 1.24); 

    model.set_window_length_rec(5);
    testTree(model, 5000, 0.9, 0.53, 2.83, 1.24); 

    model.set_window_length_rec(10);
    testTree(model, 5000, 0.9, 0.53, 2.83, 1.24); 
  }

  void testMigration() {
    Model model = Model(0);
    model.setRecombinationRate(1, false, true);
    model.set_population_number(2);
    std::vector<size_t> sample_size;
    sample_size.push_back(7);
    sample_size.push_back(3);
    model.addSampleSizes(0.0, sample_size);
    model.addSymmetricMigration(0.0, 0.5, true, true);
    model.finalize();
    testTree(model, 2000, 2.76, 1.79, 7.82, 3.86); 

    model = Model(0);
    model.setRecombinationRate(1, false, true);
    model.set_population_number(2);
    model.addSampleSizes(0.0, sample_size);
    model.addSymmetricMigration(0.0, 0.5, true, true);
    model.addSymmetricMigration(0.3, 1.1, true, true);
    model.set_window_length_seq(5);
    model.finalize();
    testTree(model, 2000, 2.24, 1.36, 6.73, 3.04); 

    model = Model(0);
    char *argv[] = { "scrm", "20", "10", "-I", "2", "10", "10", "-ma", "x", "5", "7", "x" };
    model = Param(12, argv).parse();
    testTree(model, 2000, 1.93, 1.09, 7.24, 2.56); 
  }

  void testSizeChange() {
    Model model = Model(0);
    model = Param("10 1 -r 1 100 -I 3 3 3 4 0.5 -eN 0.1 0.05 -eN 0.2 0.5").parse();
    testTree(model, 1000, 3.75, 2.68, 9.31, 5.67); 

    model.set_window_length_seq(5);
    testTree(model, 1000, 3.75, 2.68, 9.31, 5.67); 
  }


  void testGrowth() {
    Model model = Param("10 30 -r 50 50 -G 5").parse();
    testTree(model, 2500, 0.321, 0.089, 1.31, 0.28); 

    model = Param("10 30 -r 50 50 -G -0.5 -eG 0.75 2 -l 5").parse();
    testTree(model, 2500, 0.918, 0.38, 2.95, 1.00); 

    model = Param("4 30 -r 50 50 -G -2.5 -eN 1 0.25 -eG 2 0.0 -eN 2.5 0.25 -l 5").parse();
    testTree(model, 1000, 0.964, 0.34, 2.48, 0.97); 
  }

  void testSplit() {
    Model model = Model(0);
    model.setLocusLength(100);
    model.setRecombinationRate(1, false, true);
    model.set_population_number(2);
    std::vector<size_t> sample_size;
    sample_size.push_back(7);
    sample_size.push_back(3);
    model.addSampleSizes(0.0, sample_size);
    model.addSymmetricMigration(0.0, 0.5, true, true);
    model.addSingleMigrationEvent(1.0, 1, 0, 1.0, true); 
    model.addSymmetricMigration(1.0, 0.0, true, true);
    model.finalize();
    testTree(model, 2500, 1.51, 0.55, 5.20, 1.44); 

    model = Param("15 5 -r 50 50 -I 3 7 3 5 0.5 -ej 0.5 3 2 -ej 1.0 2 1").parse();
    testTree(model, 2500, 1.60, 0.54, 7.29, 1.48); 

    model.set_window_length_seq(5);
    testTree(model, 2500, 1.60, 0.54, 7.29, 1.48); 
  }

  void testMerge() {
    Model model = Param("20 10 -r 50 50 -I 2 10 10 1.5 -es 1.6 2 0.5 -eM 2.0 1").parse();
    testTree(model, 1000, 2.88, 2.26, 9.36, 4.87); 

    model.set_window_length_seq(5);
    testTree(model, 1000, 2.88, 2.26, 9.36, 4.87); 


    model = Param("20 10 -r 50 50 -I 2 10 10 1.5 -es 0.9 1 0.8 -es 1.6 2 0.5 -eM 2.0 1").parse();
    testTree(model, 1000, 3.82, 3.32, 11.39, 7.09); 

    model.set_window_length_seq(5);
    testTree(model, 1000, 3.82, 3.32, 11.39, 7.09); 
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestAlgorithm );

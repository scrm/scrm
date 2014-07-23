/*
 * Long running tests that ensure that scrm produces correct results
 * in various settings.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include <valarray>
#include <cmath>
#include "../src/forest.h"

class TestAlgorithm : public CppUnit::TestCase {

	CPPUNIT_TEST_SUITE( TestAlgorithm );

	CPPUNIT_TEST( testInitialTree );
	CPPUNIT_TEST( testARG );
	CPPUNIT_TEST( testPruning );
	CPPUNIT_TEST( testGrowth );
	CPPUNIT_TEST( testMigration );

	CPPUNIT_TEST_SUITE_END();

 private:
  MersenneTwister *rg;

  void testTree(Model &model, size_t replicates, 
                double tmrca_mean, double tmrca_sd,
                double tree_length_mean, double tree_length_sd) {
    double tmrca[6] = { 0 }, tree_length[6] = { 0 };

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
      tree_length[i] /= replicates; 
      double SE = tmrca_sd / sqrt(replicates);
      if (tmrca[i] < tmrca_mean - 3 * SE || tmrca_mean + 3 * SE < tmrca[i]) {
        std::cout << std::endl 
                  << "TMRCA outside expected range. Observed: " << tmrca[i] 
                  << " Expected: " << tmrca_mean << std::endl;  
        CPPUNIT_ASSERT( false );
      }
      SE = tree_length_sd / sqrt(replicates);
      if (tree_length[i] < tree_length_mean - 3 * SE || tree_length_mean + 3 * SE < tree_length[i]) {
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
    testTree(model, 5000, 0.8, 0.53, 2.08, 1.19); 

    model = Model(10);
    model.setRecombinationRate(0);
    testTree(model, 5000, 0.9, 0.53, 2.83, 1.24); 

    model = Model(20);
    model.setRecombinationRate(0);
    testTree(model, 5000, 0.95, 0.54, 3.55, 1.26); 
  }

  void testARG() {
    Model model = Model(5);
    model.setRecombinationRate(1, false, true);
    testTree(model, 5000, 0.8, 0.53, 2.08, 1.19); 

    model = Model(10);
    model.setRecombinationRate(1, false, true);
    testTree(model, 5000, 0.9, 0.53, 2.83, 1.24); 

    model = Model(20);
    model.setRecombinationRate(1, false, true);
    testTree(model, 5000, 0.95, 0.54, 3.55, 1.26); 
  }

  void testPruning() {
    Model model = Model(10);
    model.setRecombinationRate(1, false, true);
    model.set_exact_window_length(0);
    testTree(model, 5000, 0.9, 0.53, 2.83, 1.24); 

    model.set_exact_window_length(5);
    testTree(model, 5000, 0.9, 0.53, 2.83, 1.24); 

    model.set_exact_window_length(10);
    testTree(model, 5000, 0.9, 0.53, 2.83, 1.24); 
  }

  void testGrowth() {
    Model model = Model(10);
    model.setRecombinationRate(1, false, true);
    model.addGrowthRates(0.0, 5, true, true);
    model.finalize();
    testTree(model, 5000, 0.321, 0.089, 1.31, 0.28); 

    model.addGrowthRates(0.0, -0.5, true, true);
    model.addGrowthRates(0.75, 2, true, true);
    model.finalize();
    testTree(model, 5000, 0.918, 0.38, 2.95, 1.00); 
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
    testTree(model, 5000, 2.76, 1.79, 7.82, 3.86); 
    
    model.addSymmetricMigration(0.3, 1.1, true, true);
    model.finalize();
    std::cout << model << std::endl;
    testTree(model, 5000, 2.24, 1.36, 6.73, 3.04); 
  }

};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestAlgorithm );

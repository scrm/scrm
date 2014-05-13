/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include <valarray>
#include "../src/forest.h"

class TestAlgorithm : public CppUnit::TestCase {

	CPPUNIT_TEST_SUITE( TestAlgorithm );

	CPPUNIT_TEST( testInitialTree );
	CPPUNIT_TEST( testTreeWithGrowth );
	CPPUNIT_TEST( testTreeAfterRecombination );
	CPPUNIT_TEST( testTreeWithPruning );
	CPPUNIT_TEST( testTreeWithFullPruning );

	CPPUNIT_TEST_SUITE_END();

 private:
  MersenneTwister *rg;
  Model *model;

 public:
  void setUp() {
    model = new Model(10);
    model->set_recombination_rate(0.0001, 1000);
    rg = new MersenneTwister(78361);
  }

  void tearDown() {
    delete model;
    delete rg;
  }

  void testInitialTree() {
    double tmrca = 0;
    double tree_length = 0;
    size_t reps = 1000;

    for (size_t i = 0; i < reps; ++i) {
      Forest forest = Forest(model, rg);

      forest.buildInitialTree();
      tmrca += forest.local_root()->height() / ( 4 * model->default_pop_size );
      tree_length += forest.local_tree_length() / ( 4 * model->default_pop_size );
    }
    tmrca /= reps;        // Expectation: 0.9
    tree_length /= reps;  // Expectation: 2.83
    //std::cout << std::endl << tmrca << std::endl;
    //std::cout << tree_length << std::endl; 
    
    CPPUNIT_ASSERT( 0.85 <= tmrca && tmrca <= 0.95 );
    CPPUNIT_ASSERT( 2.75 <= tree_length && tree_length <= 2.95 );

    tmrca = 0;
    tree_length = 0;
    Model* model2 = new Model(5);

    for (size_t i = 0; i < reps; ++i) {
      Forest forest = Forest(model2, rg);

      forest.buildInitialTree();
      tmrca += forest.local_root()->height() / ( 4 * model->default_pop_size );
      tree_length += forest.local_tree_length() / ( 4 * model->default_pop_size );
    }
    tmrca /= reps;        // Expectation: 0.8
    tree_length /= reps;  // Expectation: 2.08 
    //std::cout << std::endl << tmrca << std::endl;
    //std::cout << tree_length << std::endl; 
    
    delete model2;
    CPPUNIT_ASSERT( 0.75 <= tmrca && tmrca <= 0.85 );
    CPPUNIT_ASSERT( 2.00 <= tree_length && tree_length <= 2.18 );

    tmrca = 0;
    tree_length = 0;
    Model* model3 = new Model(20);

    for (size_t i = 0; i < reps; ++i) {
      Forest forest = Forest(model3, rg);

      forest.buildInitialTree();
      tmrca += forest.local_root()->height() / ( 4 * model->default_pop_size );
      tree_length += forest.local_tree_length() / ( 4 * model->default_pop_size );
    }
    tmrca /= reps;        // Expectation: 0.95
    tree_length /= reps;  // Expectation: 3.55
    //std::cout << std::endl << tmrca << std::endl;
    //std::cout << tree_length << std::endl; 
    
    delete model3;
    CPPUNIT_ASSERT( 0.90 <= tmrca && tmrca <= 1.0 );
    CPPUNIT_ASSERT( 3.45 <= tree_length && tree_length <= 3.65 );
  }

  void testTreeAfterRecombination() {
    double tmrca[5] = { 0 };
    double tree_length[5] = { 0 };
    size_t reps = 10000;

    for (size_t i = 0; i < reps; ++i) {
      Forest forest = Forest(model, rg);
      forest.buildInitialTree();

      for (size_t j = 1; j <= 5; j++) {
        while (forest.next_base() < j*5) {
          forest.sampleNextGenealogy();
        }
        tmrca[j-1] += forest.local_root()->height() / ( 4 * model->default_pop_size );
        tree_length[j-1] += forest.local_tree_length() / ( 4 * model->default_pop_size );
      }
    }

    for (int i = 0; i < 5; ++i) {
      tmrca[i] /= reps;          // Expectation: 0.9 
      tree_length[i] /= reps;    // Expectation: 2.84 
      //std::cout << tmrca[i] << " " << tree_length[i] << std::endl; 
      CPPUNIT_ASSERT( 0.88 <= tmrca[i] && tmrca[i] <= 0.92 );
      CPPUNIT_ASSERT( 2.78 <= tree_length[i] && tree_length[i] <= 2.88 );
    }
  }

  void testTreeWithPruning() {
    double tmrca[5] = { 0 };
    double tree_length[5] = { 0 };
    size_t reps = 2000;

    model->set_exact_window_length(5);

    for (size_t i = 0; i < reps; ++i) {
      Forest forest = Forest(model, rg);
      forest.buildInitialTree();

      for (size_t j = 1; j <= 5; ++j) {
        while (forest.next_base() < j*5) {
          forest.sampleNextGenealogy();
        }
        tmrca[j-1] += forest.local_root()->height() / ( 4 * model->default_pop_size );
        tree_length[j-1] += forest.local_tree_length() / ( 4 * model->default_pop_size );
      }
    }

    for (int i = 0; i < 5; ++i) {
      tmrca[i] /= reps;          // Expectation: 0.9 
      tree_length[i] /= reps;    // Expectation: 2.84 
      //std::cout << tmrca[i] << " " << tree_length[i] << std::endl; 
      CPPUNIT_ASSERT( 0.85 <= tmrca[i] && tmrca[i] <= 0.95 );
      CPPUNIT_ASSERT( 2.74 <= tree_length[i] && tree_length[i] <= 2.94 );
    }
  }

  void testTreeWithFullPruning() {
    double tmrca[5] = { 0 };
    double tree_length[5] = { 0 };
    size_t reps = 2000;

    model->set_exact_window_length(0);

    for (size_t i = 0; i < reps; ++i) {
      Forest forest = Forest(model, rg);
      forest.buildInitialTree();

      for (size_t j = 1; j <= 5; j++) {
        while (forest.next_base() < j*5) {
          forest.sampleNextGenealogy();
        }
        tmrca[j-1] += forest.local_root()->height() / ( 4 * model->default_pop_size );
        tree_length[j-1] += forest.local_tree_length() / ( 4 * model->default_pop_size );
      }
    }

    for (int i = 0; i < 5; ++i) {
      tmrca[i] /= reps;          // Expectation: 0.9 
      tree_length[i] /= reps;    // Expectation: 2.84 
      //std::cout << tmrca[i] << " " << tree_length[i] << std::endl; 
      CPPUNIT_ASSERT( 0.85 <= tmrca[i] && tmrca[i] <= 0.95 );
      CPPUNIT_ASSERT( 2.74 <= tree_length[i] && tree_length[i] <= 2.94 );
    }
  }

  void testTreeWithGrowth() {
    double tmrca = 0 ;
    double tree_length = 0;
    size_t reps = 10000;

    model->addGrowthRates(0.0, 5, true, true);
    model->finalize();

    for (size_t i = 0; i < reps; ++i) {
      Forest forest = Forest(model, rg);
      forest.buildInitialTree();
      tmrca += forest.local_root()->height() / ( 4 * model->default_pop_size );
      tree_length += forest.local_tree_length() / ( 4 * model->default_pop_size );
    }

    tmrca /= reps;          // Expectation: 0.3209 
    tree_length /= reps;    // Expectation: 0.6727
    //std::cout << tmrca << " " << tree_length << std::endl; 
    CPPUNIT_ASSERT( 0.30 <= tmrca && tmrca <= 0.34 );
    CPPUNIT_ASSERT( 1.30 <= tree_length && tree_length <= 1.38 );
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestAlgorithm );

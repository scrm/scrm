/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include <../src/seg.h>
#include <../src/forest.h>
#include <valarray>

class TestAlgorithm : public CppUnit::TestCase {

	CPPUNIT_TEST_SUITE( TestAlgorithm );

	CPPUNIT_TEST( testInitialTree );
	CPPUNIT_TEST( testHeightChange );
	CPPUNIT_TEST( testTreeAfterRecombination );

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
    double tmrca[3] = { 0 };
    double tree_length[3] = { 0 };
    size_t reps = 2000;

    for (size_t i = 0; i < reps; ++i) {
      Forest forest = Forest(model, rg);

      int j = 0;
      forest.buildInitialTree();
      while (forest.next_base() < 5) {
        forest.sampleNextGenealogy();
      }
      tmrca[0] += forest.local_root()->height() / ( 4 * model->default_pop_size );
      tree_length[0] += forest.local_tree_length() / ( 4 * model->default_pop_size );

      while (forest.next_base() < 10) {
        forest.sampleNextGenealogy();
      }
      tmrca[1] += forest.local_root()->height() / ( 4 * model->default_pop_size );
      tree_length[1] += forest.local_tree_length() / ( 4 * model->default_pop_size );

      while (forest.next_base() < 15) {
        forest.sampleNextGenealogy();
      }
      tmrca[2] += forest.local_root()->height() / ( 4 * model->default_pop_size );
      tree_length[2] += forest.local_tree_length() / ( 4 * model->default_pop_size );
    }
    tmrca[0] /= reps;          // Expectation: 0.9 
    tmrca[1] /= reps;          // Expectation: 0.9 
    tmrca[2] /= reps;          // Expectation: 0.9 
    tree_length[0] /= reps;    // Expectation: 2.84 
    tree_length[1] /= reps;    // Expectation: 2.84 
    tree_length[2] /= reps;    // Expectation: 2.84 

    std::cout << std::endl << tmrca[0] << " " << tmrca[1] << " " << tmrca[2] << std::endl;
    std::cout << tree_length[0] <<  " " << tree_length[1] << " " << tree_length[2] << std::endl; 

    for (int i = 0; i < 3; ++i) {
      CPPUNIT_ASSERT( 0.88 <= tmrca[i] && tmrca[i] <= 0.92 );
      CPPUNIT_ASSERT( 2.80 <= tree_length[i] && tree_length[i] <= 2.90 );
    }
  }

  void testHeightChange() {
    size_t reps = 100000; 
    double lower = 0, higher = 0, equal = 0;
    double height = 0.0, new_height = 0.0;

    for (size_t i = 0; i < reps; ++i) {
      Forest forest = Forest(model, rg);
      forest.buildInitialTree();
      height = forest.local_root()->height();
      forest.sampleNextGenealogy();
      new_height = forest.local_root()->height();

      //std::cout << i << " " << height << " " << new_height << std::endl; 
      if (new_height == height) ++equal;
      if (new_height > height) ++higher;
      if (new_height < height) ++lower;
    }

    equal /= reps;
    higher /= reps; 
    lower /= reps; 

    std::cout << std::endl << equal << " " << higher << " " << lower << std::endl;
    CPPUNIT_ASSERT( 0.650 < equal && equal < 0.656 );   // E: 0.653
    CPPUNIT_ASSERT( 0.171 < lower && lower < 0.177 );   // E: 0.174
    CPPUNIT_ASSERT( 0.171 < higher && higher < 0.177 ); // E: 0.174
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestAlgorithm );

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
	CPPUNIT_TEST( testTreeAfterRecombination );

	CPPUNIT_TEST_SUITE_END();

 private:
  MersenneTwister *rg;
  Model *model;

 public:
  void setUp() {
    model = new Model(10);
    rg = new MersenneTwister(78361);
  }

  void tearDown() {
    delete model;
    delete rg;
  }


  void testInitialTree() {
    double tmrca = 0;
    double tree_length = 0;

    for (size_t i = 0; i < 1000; ++i) {
      Forest forest = Forest(model, rg);

      forest.buildInitialTree();
      tmrca += forest.local_root()->height() / ( 4 * model->default_pop_size );
      tree_length += forest.local_tree_length() / ( 4 * model->default_pop_size );
    }
    tmrca /= 1000;        // Expectation: 0.9
    tree_length /= 1000;  // Expectation: 2.84
    std::cout << std::endl << tmrca << std::endl;
    std::cout << tree_length << std::endl; 
    
    CPPUNIT_ASSERT( 0.85 <= tmrca && tmrca <= 0.95 );
    CPPUNIT_ASSERT( 2.75 <= tree_length && tree_length <= 2.95 );
  }


  void testTreeAfterRecombination() {
    double tmrca = 0;
    double tree_length = 0;

    for (size_t i = 0; i < 1000; ++i) {
      model->set_recombination_rate(10, 100);
      Forest forest = Forest(model, rg);

      forest.buildInitialTree();
      for (int j=0; j < 20; ++j) forest.sampleNextGenealogy();
      tmrca += forest.local_root()->height() / ( 4 * model->default_pop_size );
      tree_length += forest.local_tree_length() / ( 4 * model->default_pop_size );
    }
    tmrca /= 1000;          // Expectation: 0.9 
    tree_length /= 1000;    // Expectation: 2.84 

    std::cout << std::endl << tmrca << std::endl;
    std::cout << tree_length << std::endl; 
    CPPUNIT_ASSERT( 0.85 <= tmrca && tmrca <= 0.95 );
    CPPUNIT_ASSERT( 2.75 <= tree_length && tree_length <= 2.95 );
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestAlgorithm );

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

 public:
  void setUp() {
  }

  void tearDown() {
  }

  void testInitialTree() {
    double tmrca = 0;
    double tree_length = 0;

    for (size_t i = 0; i < 1000; ++i) {
      Model* model = new Model(10);
      rg = new MersenneTwister(i);
      Forest forest = Forest(model, rg);

      forest.buildInitialTree();
      tmrca += forest.local_root()->height() / ( 2 * model->default_pop_size );
      tree_length += forest.local_tree_length() / ( 2 * model->default_pop_size );
    }
    tmrca /= 1000; // Expectation: 1.8
    tree_length /= 1000; // Expectation: 5.66
    //std::cout << tmrca << std::endl;
    //std::cout << tree_length << std::endl; 
    
    CPPUNIT_ASSERT( 1.75 <= tmrca && tmrca <= 1.85 );
    CPPUNIT_ASSERT( 5.5 <= tree_length && tree_length <= 5.8 );

  }

  void testTreeAfterRecombination() {
    double tmrca = 0;
    double tree_length = 0;

    for (size_t i = 0; i < 1000; ++i) {
      Model* model = new Model(10);
      model->set_recombination_rate(10, 100);
      rg = new MersenneTwister(i+10000);
      Forest forest = Forest(model, rg);

      forest.buildInitialTree();
      for (int j=0; j < 20; ++j) forest.sampleNextGenealogy();
      tmrca += forest.local_root()->height() / ( 2 * model->default_pop_size );
      tree_length += forest.local_tree_length() / ( 2 * model->default_pop_size );
    }
    tmrca /= 1000;          // Expectation: 1.8
    tree_length /= 1000;    // Expectation: 5.66

    std::cout << std::endl << tmrca << std::endl;
    std::cout << tree_length << std::endl; 
    CPPUNIT_ASSERT( 1.75 <= tmrca && tmrca <= 1.85 );
    CPPUNIT_ASSERT( 5.5 <= tree_length && tree_length <= 5.8 );
  }

};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestAlgorithm );

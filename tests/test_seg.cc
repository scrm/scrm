/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include <../src/seg.h>
#include <../src/forest.h>
#include <valarray>

class TestSeg : public CppUnit::TestCase {

	CPPUNIT_TEST_SUITE( TestSeg );

	CPPUNIT_TEST( testTraversal );

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
    delete forest->writable_model();
    delete forest;
    delete rg;
  }

  void testTraversal() {
    std::valarray <int>haplotype(4);
    traversal(forest->nodes()->at(4), haplotype);
    CPPUNIT_ASSERT_EQUAL( 1, haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( 1, haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( 0, haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( 0, haplotype[3] );

    haplotype *= 0;
    traversal(forest->nodes()->at(5), haplotype);
    CPPUNIT_ASSERT_EQUAL( 0, haplotype[0] );
    CPPUNIT_ASSERT_EQUAL( 0, haplotype[1] );
    CPPUNIT_ASSERT_EQUAL( 1, haplotype[2] );
    CPPUNIT_ASSERT_EQUAL( 1, haplotype[3] );
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestSeg );

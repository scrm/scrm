#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../../src/contemporaries_container.h"
#include "../../src/random/mersenne_twister.h"

class TestContemporariesContainer : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestContemporariesContainer );

  CPPUNIT_TEST( add );
  CPPUNIT_TEST( remove );
  CPPUNIT_TEST( clear );
  CPPUNIT_TEST( iterator );
  CPPUNIT_TEST( sample );

  CPPUNIT_TEST_SUITE_END();

 private:
  MersenneTwister *rg;
  Node *node1, *node2, *node3, *node4;

 public:
  void setUp() {
    rg = new MersenneTwister(5);
    node1 = new Node(5);
    node2 = new Node(10);
    node2->set_population(1);
    node3 = new Node(10);
    node3->set_population(2);
    node4 = new Node(15);
    node4->set_population(2);
  }

  void tearDown() {
    delete rg;
    delete node1, node2, node3, node4;
  }

  void add() {
    ContemporariesContainer cc = ContemporariesContainer(3, 10, rg);
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(2) );

    cc.add(node1);
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(2) );

    cc.add(node2);
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(2) );

    cc.add(node3);
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(2) );

    cc.add(node4);
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, cc.size(2) );
  }

  void remove() {
    // Test Remove
    ContemporariesContainer cc = ContemporariesContainer(3, 10, rg);
    cc.add(node1);
    cc.add(node2);
    cc.add(node3);
    cc.add(node4);
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, cc.size(2) );

    cc.remove(node1);
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)2, cc.size(2) );

    cc.remove(node3);
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(2) );

    cc.remove(node4);
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)1, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(2) );

    cc.remove(node2);
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(2) );
  }

  void clear() {
    ContemporariesContainer cc = ContemporariesContainer(3, 10, rg);
    cc.add(node1);
    cc.add(node2);
    cc.add(node3);
    cc.add(node4);

    cc.clear();
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(2) );
  }

  void iterator() {
    ContemporariesContainer cc = ContemporariesContainer(3, 10, rg);
    cc.add(node1);
    cc.add(node2);
    cc.add(node3);
    cc.add(node4);

    size_t count = 0;
    for (auto it = cc.begin(0); it != cc.end(0); ++it) {
      CPPUNIT_ASSERT( *it == node1 );
      ++count;
    }
    CPPUNIT_ASSERT_EQUAL((size_t)1, count);

    for (auto it = cc.begin(2); it != cc.end(2); ++it) {
      CPPUNIT_ASSERT( *it == node4 || *it == node3 );
      ++count;
    }
    CPPUNIT_ASSERT_EQUAL((size_t)3, count);

    for (auto it = cc.begin(1); it != cc.end(1); ++it) {
      CPPUNIT_ASSERT( *it == node2 );
      ++count;
    }
    CPPUNIT_ASSERT_EQUAL((size_t)4, count);
  }

  void sample() {
    ContemporariesContainer cc = ContemporariesContainer(3, 10, rg);
    cc.add(node1);
    cc.add(node2);
    cc.add(node3);
    cc.add(node4);

    for (size_t i = 0; i < 1000; ++i) {
      CPPUNIT_ASSERT_EQUAL(node1, cc.sample(0));
      CPPUNIT_ASSERT_EQUAL(node2, cc.sample(1));
    } 
  
    double count = 0;
    for (size_t i = 0; i < 10000; ++i) {
      Node* node = cc.sample(2);
      CPPUNIT_ASSERT( node == node3 || node == node4 );
      if (node == node3) ++count;
    } 
    count /= 10000;
    CPPUNIT_ASSERT( 0.49 < count && count < 0.51 ); 
  }
};
CPPUNIT_TEST_SUITE_REGISTRATION( TestContemporariesContainer );

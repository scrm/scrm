#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>

#include "../../src/contemporaries_container.h"
#include "../../src/node_container.h"
#include "../../src/random/mersenne_twister.h"

class TestContemporariesContainer : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestContemporariesContainer );

  CPPUNIT_TEST( add );
  CPPUNIT_TEST( remove );
  CPPUNIT_TEST( clear );
  CPPUNIT_TEST( iterator );
  CPPUNIT_TEST( sample );
  CPPUNIT_TEST( buffer );
  CPPUNIT_TEST( empty );

  CPPUNIT_TEST_SUITE_END();

 private:
  MersenneTwister *rg;
  NodeContainer *nc;
  Node *node1, *node2, *node3, *node4;

 public:
  void setUp() {
    rg = new MersenneTwister(5);
    nc = new NodeContainer();
    node1 = nc->createNode(5);
    node2 = nc->createNode(10);
    node2->set_population(1);
    node3 = nc->createNode(10);
    node3->set_population(2);
    node4 = nc->createNode(15);
    node4->set_population(2);
  }

  void tearDown() {
    delete rg, nc;
  }

  void add() { 
    // Vector
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

    // Set
    cc = ContemporariesContainer(3, 1000, rg);
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
    // Vector
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

    // Set
    cc = ContemporariesContainer(3, 1000, rg);
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
    
    cc = ContemporariesContainer(3, 1000, rg);
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
    // Vector
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

    // Set
    cc = ContemporariesContainer(3, 1000, rg);
    cc.add(node1);
    cc.add(node2);
    cc.add(node3);
    cc.add(node4);

    count = 0;
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
    // Vector
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

    // Set
    cc = ContemporariesContainer(3, 1000, rg);
    cc.add(node1);
    cc.add(node2);
    cc.add(node3);
    cc.add(node4);

    for (size_t i = 0; i < 1000; ++i) {
      CPPUNIT_ASSERT_EQUAL(node1, cc.sample(0));
      CPPUNIT_ASSERT_EQUAL(node2, cc.sample(1));
    } 

    count = 0;
    for (size_t i = 0; i < 10000; ++i) {
      Node* node = cc.sample(2);
      CPPUNIT_ASSERT( node == node3 || node == node4 );
      if (node == node3) ++count;
    } 
    count /= 10000;
    CPPUNIT_ASSERT( 0.49 < count && count < 0.51 ); 
  }

  void buffer() {
    // Vector
    ContemporariesContainer cc = ContemporariesContainer(3, 10, rg);
    cc.add(node1);
    cc.add(node2);
    cc.add(node3);
    cc.add(node4);

    cc.buffer(17.5);
    CPPUNIT_ASSERT_EQUAL( 17.5, cc.buffer_time() );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(2) );

    size_t count = 0;
    for (auto it = cc.buffer_begin(0); it != cc.buffer_end(0); ++it) {
      CPPUNIT_ASSERT( *it == node1 );
      ++count;
    }
    CPPUNIT_ASSERT_EQUAL((size_t)1, count);

    for (auto it = cc.buffer_begin(2); it != cc.buffer_end(2); ++it) {
      CPPUNIT_ASSERT( *it == node4 || *it == node3 );
      ++count;
    }
    CPPUNIT_ASSERT_EQUAL((size_t)3, count);

    for (auto it = cc.buffer_begin(1); it != cc.buffer_end(1); ++it) {
      CPPUNIT_ASSERT( *it == node2 );
      ++count;
    }
    CPPUNIT_ASSERT_EQUAL((size_t)4, count);

    cc.buffer(20.2);
    CPPUNIT_ASSERT_EQUAL( 20.2, cc.buffer_time() );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(2) );
    count = 0;
    for (auto it = cc.buffer_begin(0); it != cc.buffer_end(0); ++it) ++count;
    for (auto it = cc.buffer_begin(1); it != cc.buffer_end(1); ++it) ++count;
    for (auto it = cc.buffer_begin(2); it != cc.buffer_end(2); ++it) ++count;
    CPPUNIT_ASSERT_EQUAL((size_t)0, count);

    // Set
    cc = ContemporariesContainer(3, 1000, rg);
    cc.add(node1);
    cc.add(node2);
    cc.add(node3);
    cc.add(node4);

    cc.buffer(17.5);
    CPPUNIT_ASSERT_EQUAL( 17.5, cc.buffer_time() );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(2) );

    count = 0;
    for (auto it = cc.buffer_begin(0); it != cc.buffer_end(0); ++it) {
      CPPUNIT_ASSERT( *it == node1 );
      ++count;
    }
    CPPUNIT_ASSERT_EQUAL((size_t)1, count);

    for (auto it = cc.buffer_begin(2); it != cc.buffer_end(2); ++it) {
      CPPUNIT_ASSERT( *it == node4 || *it == node3 );
      ++count;
    }
    CPPUNIT_ASSERT_EQUAL((size_t)3, count);

    for (auto it = cc.buffer_begin(1); it != cc.buffer_end(1); ++it) {
      CPPUNIT_ASSERT( *it == node2 );
      ++count;
    }
    CPPUNIT_ASSERT_EQUAL((size_t)4, count);

    cc.buffer(20.2);
    CPPUNIT_ASSERT_EQUAL( 20.2, cc.buffer_time() );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(0) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(1) );
    CPPUNIT_ASSERT_EQUAL( (size_t)0, cc.size(2) );
    count = 0;
    for (auto it = cc.buffer_begin(0); it != cc.buffer_end(0); ++it) ++count;
    for (auto it = cc.buffer_begin(1); it != cc.buffer_end(1); ++it) ++count;
    for (auto it = cc.buffer_begin(2); it != cc.buffer_end(2); ++it) ++count;
    CPPUNIT_ASSERT_EQUAL((size_t)0, count);
  }

  void empty() {
    // Vector 
    ContemporariesContainer cc = ContemporariesContainer(3, 10, rg);
    CPPUNIT_ASSERT( cc.empty() );
    cc.add(node1);
    CPPUNIT_ASSERT( !cc.empty() );

    cc.clear();
    CPPUNIT_ASSERT( cc.empty() );
    cc.add(node2);
    CPPUNIT_ASSERT( !cc.empty() );

    cc.clear();
    CPPUNIT_ASSERT( cc.empty() );
    cc.add(node3);
    CPPUNIT_ASSERT( !cc.empty() );

    // Set
    cc = ContemporariesContainer(3, 1000, rg);
    CPPUNIT_ASSERT( cc.empty() );
    cc.add(node1);
    CPPUNIT_ASSERT( !cc.empty() );

    cc.clear();
    CPPUNIT_ASSERT( cc.empty() );
    cc.add(node2);
    CPPUNIT_ASSERT( !cc.empty() );

    cc.clear();
    CPPUNIT_ASSERT( cc.empty() );
    cc.add(node3);
    CPPUNIT_ASSERT( !cc.empty() );
  }
};
CPPUNIT_TEST_SUITE_REGISTRATION( TestContemporariesContainer );

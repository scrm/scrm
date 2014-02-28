/*
 * Unit tests for NodeContainer
 * If the container passes this tests, it should also do for the remaining code
 * */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/forest.h"


class TestNodeContainer : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestNodeContainer );

  CPPUNIT_TEST( testGet );
  CPPUNIT_TEST( testAdd );
  CPPUNIT_TEST( testNodeIterator );
  CPPUNIT_TEST( testNodeIteratorHeight );
  CPPUNIT_TEST( testReverseIterator );
  CPPUNIT_TEST( testRemove );
  CPPUNIT_TEST( testMove );
  CPPUNIT_TEST( testCopyConstructor );

  CPPUNIT_TEST_SUITE_END();

 private:
    Node *node1;
    Node *node2;
    Node *node3;
    NodeContainer nc;

 public:
  void setUp() {
    node1 = new Node(1, 1);
    node2 = new Node(2, 2);
    node3 = new Node(3, 0);
    nc = NodeContainer();
    nc.add(node1);
    nc.add(node2);
    nc.add(node3);
  }

  void tearDown() {
    nc.clear();
  }

  void testGet() {
    CPPUNIT_ASSERT( nc.get(0) == node1 );
    CPPUNIT_ASSERT( nc.get(1) == node2 );
    CPPUNIT_ASSERT( nc.get(2) == node3 );
  }

  void testAdd() {
    CPPUNIT_ASSERT( nc.get(0) == node1 );
    CPPUNIT_ASSERT( nc.get(1) == node2 );
    CPPUNIT_ASSERT( nc.get(2) == node3 );
    CPPUNIT_ASSERT( nc.first() == node1 );
    CPPUNIT_ASSERT( nc.last() == node3 );
    CPPUNIT_ASSERT_EQUAL( (size_t)3, nc.size() );

    // Add new first node
    Node* node = new Node(0);
    nc.add(node);
    CPPUNIT_ASSERT_EQUAL( (size_t)4, nc.size() );
    CPPUNIT_ASSERT( nc.get(0) == node );
    CPPUNIT_ASSERT( nc.first() == node );

    // Add new last node
    node = new Node(10);
    nc.add(node);
    CPPUNIT_ASSERT_EQUAL( (size_t)5, nc.size() );
    CPPUNIT_ASSERT( nc.get(4) == node );
    CPPUNIT_ASSERT( nc.last() == node );

    // Add something in between 
    node = new Node(3);
    nc.add(node);
    CPPUNIT_ASSERT_EQUAL( (size_t)6, nc.size() );
    CPPUNIT_ASSERT( nc.get(4) == node );
  }

  void testRemove() {
    Node *node = new Node(2.5);
    nc.add(node);
    nc.remove(node);
    CPPUNIT_ASSERT( nc.size() == 3 );
    CPPUNIT_ASSERT( nc.get(0) == node1 );
    CPPUNIT_ASSERT( nc.get(1) == node2 );
    CPPUNIT_ASSERT( nc.get(2) == node3 );

    nc.remove(node3);
    CPPUNIT_ASSERT( nc.size() == 2 );
    CPPUNIT_ASSERT( nc.last() == node2 );

    nc.remove(node1);
    CPPUNIT_ASSERT( nc.size() == 1 );
    CPPUNIT_ASSERT( nc.first() == node2 );

    nc.remove(node2);
    CPPUNIT_ASSERT( nc.size() == 0 );
    CPPUNIT_ASSERT( nc.first() == NULL );
    CPPUNIT_ASSERT( nc.last() == NULL );
  }


  void testMove() {
    Node* node4 = new Node(4);
    nc.add(node4);

    // Middle -> Middle (forwards)
    nc.move(node2, 3.5);
    CPPUNIT_ASSERT( nc.get(0) == node1 );
    CPPUNIT_ASSERT( nc.get(1) == node3 );
    CPPUNIT_ASSERT( nc.get(2) == node2 );
    CPPUNIT_ASSERT( nc.get(3) == node4 );
    CPPUNIT_ASSERT( node2->height() == 3.5 );

    // Middle -> Middle (backwards)
    nc.move(node2, 2);
    CPPUNIT_ASSERT( nc.get(0) == node1 );
    CPPUNIT_ASSERT( nc.get(1) == node2 );
    CPPUNIT_ASSERT( nc.get(2) == node3 );
    CPPUNIT_ASSERT( nc.get(3) == node4 );
    CPPUNIT_ASSERT( node2->height() == 2 );

    // start -> end
    nc.move(node1, 10);
    CPPUNIT_ASSERT( nc.get(0) == node2 );
    CPPUNIT_ASSERT( nc.get(1) == node3 );
    CPPUNIT_ASSERT( nc.get(2) == node4 );
    CPPUNIT_ASSERT( nc.get(3) == node1 );
    CPPUNIT_ASSERT( nc.first() == node2 );
    CPPUNIT_ASSERT( nc.last() == node1 );
    CPPUNIT_ASSERT( node1->height() == 10 );
        
    // end -> start
    nc.move(node1, 1);
    CPPUNIT_ASSERT( nc.get(0) == node1 );
    CPPUNIT_ASSERT( nc.get(1) == node2 );
    CPPUNIT_ASSERT( nc.get(2) == node3 );
    CPPUNIT_ASSERT( nc.get(3) == node4 );
    CPPUNIT_ASSERT( nc.first() == node1 );
    CPPUNIT_ASSERT( nc.last() == node4 );
    CPPUNIT_ASSERT( node1->height() == 1 );

    // start -> start
    nc.move(node1, 0.5);
    CPPUNIT_ASSERT( nc.get(0) == node1 );
    CPPUNIT_ASSERT( nc.get(1) == node2 );
    CPPUNIT_ASSERT( nc.get(2) == node3 );
    CPPUNIT_ASSERT( nc.get(3) == node4 );
    CPPUNIT_ASSERT( nc.first() == node1 );
    CPPUNIT_ASSERT( nc.last() == node4 );
    CPPUNIT_ASSERT( node1->height() == 0.5 );

    // end -> end
    nc.move(node4, 11);
    CPPUNIT_ASSERT( nc.get(0) == node1 );
    CPPUNIT_ASSERT( nc.get(1) == node2 );
    CPPUNIT_ASSERT( nc.get(2) == node3 );
    CPPUNIT_ASSERT( nc.get(3) == node4 );
    CPPUNIT_ASSERT( nc.first() == node1 );
    CPPUNIT_ASSERT( nc.last() == node4 );
    CPPUNIT_ASSERT( node4->height() == 11 );
  }


  void testNodeIterator() {
    NodeIterator it = nc.iterator();
    CPPUNIT_ASSERT( *it == node1 );
    CPPUNIT_ASSERT( it.good() );
    CPPUNIT_ASSERT( it.good() && it++ == node1 );
    CPPUNIT_ASSERT( it.good() && it++ == node2 );
    CPPUNIT_ASSERT( it.good() && it++ == node3 );
    CPPUNIT_ASSERT( !it.good() );
    CPPUNIT_ASSERT_THROW( ++it, std::out_of_range );

    it = nc.iterator();
    int i = 0;
    for (it = nc.iterator(); it.good(); ++it) {
      (*it)->make_local();
      ++i;
    }
    CPPUNIT_ASSERT( i == 3 );

    it = nc.iterator();
    ++it; --it;
    CPPUNIT_ASSERT( it.good() && *it == node1 );
    ++it; ++it; --it;
    CPPUNIT_ASSERT( it.good() && *it == node2 );
    --it; --it;
    CPPUNIT_ASSERT_THROW( --it, std::out_of_range );
  }
  

  void testReverseIterator() {
    ReverseConstNodeIterator it = nc.reverse_iterator();
    CPPUNIT_ASSERT( *it == node3 );
    CPPUNIT_ASSERT( it.good() );
    CPPUNIT_ASSERT( it.good() && it++ == node3 );
    CPPUNIT_ASSERT( it.good() && it++ == node2 );
    CPPUNIT_ASSERT( it.good() && it++ == node1 );
    CPPUNIT_ASSERT( !it.good() );
    CPPUNIT_ASSERT_THROW( ++it, std::out_of_range );
  }
  

  void testNodeIteratorHeight() {
    NodeIterator it = nc.iterator();
    CPPUNIT_ASSERT( it.height() == 1 );
    ++it; CPPUNIT_ASSERT( it.height() == 2 );
    ++it; CPPUNIT_ASSERT( it.height() == 3 );
    ++it; CPPUNIT_ASSERT( it.height() == DBL_MAX );
  }

  void testCopyConstructor() {
    NodeContainer nc2 = NodeContainer(nc);
    CPPUNIT_ASSERT_EQUAL( (size_t)3, nc2.size() );
    CPPUNIT_ASSERT( nc2.sorted() );
    CPPUNIT_ASSERT( nc2.unsorted_node_ == NULL );
    for (auto it = nc.iterator(); it.good(); ++it) {
      CPPUNIT_ASSERT( (*it)->label() <= 2 );
    }
  }
};

//Uncomment this to make_local the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestNodeContainer );

/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/forest.h"


class TestNodeContainer : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestNodeContainer );

  CPPUNIT_TEST( testGet );
  CPPUNIT_TEST( testAdd );
  CPPUNIT_TEST( testSize );
  CPPUNIT_TEST( testRemove );
  CPPUNIT_TEST( testNodeIterator );

  CPPUNIT_TEST_SUITE_END();

 private:

 public:
  void setUp() {}
  void tearDown() {}

  void testGet() {
    NodeContainer nc = NodeContainer();
    Node *node1 = new Node(1);
    nc.add(node1);
    CPPUNIT_ASSERT( nc.get(0) == node1 );
  }

  void testAdd() {
    Node *node1 = new Node(1);
    Node *node2 = new Node(2);
    Node *node3 = new Node(3);
    NodeContainer nc = NodeContainer();

    //Sorted
    nc.add(node2);
    nc.add(node3);
    nc.add(node1);
    CPPUNIT_ASSERT( nc.get(0) == node1 );
    CPPUNIT_ASSERT( nc.get(1) == node2 );
    CPPUNIT_ASSERT( nc.get(2) == node3 );
    
    /*Unsorted
    node1 = new Node(1);
    node2 = new Node(2);
    nc.add(node2, false);
    nc.add(node1, false);
    CPPUNIT_ASSERT( nc.get(0, false) == node2 );
    CPPUNIT_ASSERT( nc.get(1, false) == node1 );
    */
  }

  void testSize() {
    NodeContainer nc = NodeContainer();
    nc.add(new Node(1));
    CPPUNIT_ASSERT( nc.size() == 1 );
    nc.add(new Node(1));
    CPPUNIT_ASSERT( nc.size() == 2 );
    nc.add(new Node(1), false);
    CPPUNIT_ASSERT( nc.size() == 3 );
  }

  void testRemove() {
    NodeContainer nc = NodeContainer();
    nc.add(new Node(1));
    nc.add(new Node(2));

    Node *node = new Node(3);
    nc.add(node);
    nc.remove(node);
    CPPUNIT_ASSERT( nc.size() == 2 );

    /*Node *node2 = new Node(3);
    nc.add(node2, false);
    nc.remove(node2);
    CPPUNIT_ASSERT( nc.size() == 2 );
    */
  }

  void testNodeIterator() {
    Node *node1 = new Node(1);
    Node *node2 = new Node(2);
    Node *node3 = new Node(3);
    NodeContainer nc = NodeContainer();
    nc.add(node1);
    nc.add(node2);
    nc.add(node3);

    //Create the ultimate root
    Node* ultimate_root = new Node(FLT_MAX);
    ultimate_root->deactivate();
    nc.add(ultimate_root);

    NodeIterator it = nc.interator();
    CPPUNIT_ASSERT( it++ == node1 );
    CPPUNIT_ASSERT( it++ == node2 );
    CPPUNIT_ASSERT( it++ == node3 );

    it = nc.interator();
    int i = 0;
    while (it.good()) {
      ++i;
      ++it;
    }
    CPPUNIT_ASSERT( i == 3 );

    for (it = nc.interator(); it.good(); ++it) {
      ++i;
    }
    CPPUNIT_ASSERT( i == 6 );
  }
};

//Uncomment this to activate the test
CPPUNIT_TEST_SUITE_REGISTRATION( TestNodeContainer );

/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/node.h"
#include "../src/forest.h"

class TestNode : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestNode );

  CPPUNIT_TEST( testGettersAndSetters );
  CPPUNIT_TEST( testIsRoot );
  CPPUNIT_TEST( testInSample );
  CPPUNIT_TEST( testSamplesBelow );
  CPPUNIT_TEST( testLengthBelow );
  CPPUNIT_TEST( testNumberOfChildren );

  CPPUNIT_TEST_SUITE_END();

 private:
  Forest *forest;
  ConstantGenerator *rg;

 public:
  void setUp() {
    rg = new ConstantGenerator();
    forest = new Forest(new Model(0), rg);
    forest->createExampleTree();
  }

  void tearDown() {
    delete forest;
  }

  void testGettersAndSetters() {
    Node node1, node2;

    //height
    node1.set_height(1);
    CPPUNIT_ASSERT( node1.height() == 1 );

    //parent
    node2.set_parent(&node1);
    CPPUNIT_ASSERT( node2.parent()->height() == 1 );

    //Children
    node2.set_second_child(&node1);
    CPPUNIT_ASSERT( node2.second_child()->height() == 1 );
    node2.set_first_child(&node1);
    CPPUNIT_ASSERT( node2.first_child()->height() == 1 );

    //local
    node1.set_local(true);
    node2.set_local(false);
    CPPUNIT_ASSERT( node1.local() && !node2.local() );
    node1.make_nonlocal(1);
    node2.make_local();
    CPPUNIT_ASSERT( (!node1.local()) && node2.local() );
  }

  void testIsRoot(){
    CPPUNIT_ASSERT( !forest->nodes()->get(0)->is_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(1)->is_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(2)->is_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(3)->is_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(4)->is_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(5)->is_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(6)->is_root() );
    CPPUNIT_ASSERT( forest->nodes()->get(7)->is_root() );
    CPPUNIT_ASSERT( forest->nodes()->get(8)->is_root() );
  }

  void testInSample(){
    CPPUNIT_ASSERT( forest->nodes()->get(0)->in_sample() );
    CPPUNIT_ASSERT( forest->nodes()->get(1)->in_sample() );
    CPPUNIT_ASSERT( forest->nodes()->get(2)->in_sample() );
    CPPUNIT_ASSERT( forest->nodes()->get(3)->in_sample() );
    CPPUNIT_ASSERT( !forest->nodes()->get(4)->in_sample() );
    CPPUNIT_ASSERT( !forest->nodes()->get(5)->in_sample() );
    CPPUNIT_ASSERT( !forest->nodes()->get(6)->in_sample() );
    CPPUNIT_ASSERT( !forest->nodes()->get(7)->in_sample() );
    CPPUNIT_ASSERT( !forest->nodes()->get(8)->in_sample() );
  }

  void testSamplesBelow(){
    CPPUNIT_ASSERT( forest->nodes()->get(0)->samples_below() == 1 );
    CPPUNIT_ASSERT( forest->nodes()->get(1)->samples_below() == 1 );
    CPPUNIT_ASSERT( forest->nodes()->get(2)->samples_below() == 1 );
    CPPUNIT_ASSERT( forest->nodes()->get(3)->samples_below() == 1 );
    CPPUNIT_ASSERT( forest->nodes()->get(4)->samples_below() == 2 );
    CPPUNIT_ASSERT( forest->nodes()->get(5)->samples_below() == 2 );
    CPPUNIT_ASSERT( forest->nodes()->get(6)->samples_below() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(7)->samples_below() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(8)->samples_below() == 4 );
  }

  void testLengthBelow(){
    CPPUNIT_ASSERT( forest->nodes()->get(0)->length_below() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(1)->length_below() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(2)->length_below() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(3)->length_below() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(4)->length_below() == 2 );
    CPPUNIT_ASSERT( forest->nodes()->get(5)->length_below() == 6 );
    CPPUNIT_ASSERT( forest->nodes()->get(6)->length_below() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(7)->length_below() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(8)->length_below() == 24 );
  }

  void testNumberOfChildren(){
    CPPUNIT_ASSERT( forest->nodes()->get(0)->numberOfChildren() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(1)->numberOfChildren() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(2)->numberOfChildren() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(3)->numberOfChildren() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(4)->numberOfChildren() == 2 );
    CPPUNIT_ASSERT( forest->nodes()->get(5)->numberOfChildren() == 2 );
    CPPUNIT_ASSERT( forest->nodes()->get(6)->numberOfChildren() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(7)->numberOfChildren() == 1 );
    CPPUNIT_ASSERT( forest->nodes()->get(8)->numberOfChildren() == 2 );
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestNode );

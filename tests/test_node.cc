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
  CPPUNIT_TEST( testIsFake );
  CPPUNIT_TEST( testIsUltimateRoot );
  CPPUNIT_TEST( testIsRoot );
  CPPUNIT_TEST( testInSample );
  CPPUNIT_TEST( testSamplesBelow );
  CPPUNIT_TEST( testLengthBelow );

  CPPUNIT_TEST_SUITE_END();

 private:
  Forest *forest;
  ConstantGenerator *rg;

 public:
  void setUp() {
    rg = new ConstantGenerator();
    forest = new Forest(Model(0), rg);
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
    node2.set_higher_child(&node1);
    CPPUNIT_ASSERT( node2.higher_child()->height() == 1 );
    node2.set_lower_child(&node1);
    CPPUNIT_ASSERT( node2.lower_child()->height() == 1 );

    //active
    node1.set_active(true);
    node2.set_active(false);
    CPPUNIT_ASSERT( node1.active() && !node2.active() );
    node1.deactivate(1);
    node2.activate();
    CPPUNIT_ASSERT( (!node1.active()) && node2.active() );
  }

  void testIsRoot(){
    CPPUNIT_ASSERT( !forest->nodes()->get(0)->is_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(1)->is_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(2)->is_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(3)->is_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(4)->is_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(5)->is_root() );
    CPPUNIT_ASSERT( forest->nodes()->get(6)->is_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(7)->is_root() );
  }

  void testIsUltimateRoot(){
    CPPUNIT_ASSERT( !forest->nodes()->get(0)->is_ultimate_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(1)->is_ultimate_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(2)->is_ultimate_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(3)->is_ultimate_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(4)->is_ultimate_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(5)->is_ultimate_root() );
    CPPUNIT_ASSERT( !forest->nodes()->get(6)->is_ultimate_root() );
    CPPUNIT_ASSERT( forest->nodes()->get(7)->is_ultimate_root() );
  }

  void testIsFake(){
    CPPUNIT_ASSERT( !forest->nodes()->get(0)->is_fake() );
    CPPUNIT_ASSERT( !forest->nodes()->get(1)->is_fake() );
    CPPUNIT_ASSERT( !forest->nodes()->get(2)->is_fake() );
    CPPUNIT_ASSERT( !forest->nodes()->get(3)->is_fake() );
    CPPUNIT_ASSERT( !forest->nodes()->get(4)->is_fake() );
    CPPUNIT_ASSERT( !forest->nodes()->get(5)->is_fake() );
    CPPUNIT_ASSERT( !forest->nodes()->get(6)->is_fake() );
    CPPUNIT_ASSERT( forest->nodes()->get(7)->is_fake() );
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
  }

  void testSamplesBelow(){
    CPPUNIT_ASSERT( forest->nodes()->get(0)->samples_below() == 1 );
    CPPUNIT_ASSERT( forest->nodes()->get(1)->samples_below() == 1 );
    CPPUNIT_ASSERT( forest->nodes()->get(2)->samples_below() == 1 );
    CPPUNIT_ASSERT( forest->nodes()->get(3)->samples_below() == 1 );
    CPPUNIT_ASSERT( forest->nodes()->get(4)->samples_below() == 2 );
    CPPUNIT_ASSERT( forest->nodes()->get(5)->samples_below() == 2 );
    CPPUNIT_ASSERT( forest->nodes()->get(6)->samples_below() == 4 );
  }

  void testLengthBelow(){
    CPPUNIT_ASSERT( forest->nodes()->get(0)->length_below() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(1)->length_below() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(2)->length_below() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(3)->length_below() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(4)->length_below() == 2 );
    CPPUNIT_ASSERT( forest->nodes()->get(5)->length_below() == 6 );
    CPPUNIT_ASSERT( forest->nodes()->get(6)->length_below() == 24 );
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( TestNode );

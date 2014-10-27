/*
 * A sample test case which can be used as a template.
 */
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../../src/node.h"
#include "../../src/forest.h"
#include "../../src/random/constant_generator.h"

class TestNode : public CppUnit::TestCase {

  CPPUNIT_TEST_SUITE( TestNode );

  CPPUNIT_TEST( testGettersAndSetters );
  CPPUNIT_TEST( testIsRoot );
  CPPUNIT_TEST( testInSample );
  CPPUNIT_TEST( testSamplesBelow );
  CPPUNIT_TEST( testLengthBelow );
  CPPUNIT_TEST( testCountChildren );
  CPPUNIT_TEST( testLocalNavigation );

  CPPUNIT_TEST_SUITE_END();

 private:
  Forest *forest;
  ConstantGenerator *rg;
  Model *model;

 public:
  void setUp() {
    rg = new ConstantGenerator();
    model = new Model(0);
    forest = new Forest(model, rg);
    forest->createExampleTree();
  }

  void tearDown() {
    delete forest;
    delete rg;
    delete model;
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
    node1.make_local();
    node2.make_nonlocal(1);
    CPPUNIT_ASSERT( node1.local() && !node2.local() );
    node1.make_nonlocal(1);
    node2.make_local();
    CPPUNIT_ASSERT( (!node1.local()) && node2.local() );

    //population
    CPPUNIT_ASSERT_EQUAL( (size_t)0, node1.population() );
    node1.set_population(1);
    CPPUNIT_ASSERT_EQUAL( (size_t)1, node1.population() );
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

  void testCountChildren(){
    CPPUNIT_ASSERT( forest->nodes()->get(0)->countChildren() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(1)->countChildren() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(2)->countChildren() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(3)->countChildren() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(4)->countChildren() == 2 );
    CPPUNIT_ASSERT( forest->nodes()->get(5)->countChildren() == 2 );
    CPPUNIT_ASSERT( forest->nodes()->get(6)->countChildren() == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(7)->countChildren() == 1 );
    CPPUNIT_ASSERT( forest->nodes()->get(8)->countChildren() == 2 );

    CPPUNIT_ASSERT( forest->nodes()->get(0)->countChildren(true) == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(1)->countChildren(true) == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(2)->countChildren(true) == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(3)->countChildren(true) == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(4)->countChildren(true) == 2 );
    CPPUNIT_ASSERT( forest->nodes()->get(5)->countChildren(true) == 2 );
    CPPUNIT_ASSERT( forest->nodes()->get(6)->countChildren(true) == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(7)->countChildren(true) == 0 );
    CPPUNIT_ASSERT( forest->nodes()->get(8)->countChildren(true) == 2 );

    forest->nodes()->at(4)->make_nonlocal(1.0);
    CPPUNIT_ASSERT( forest->nodes()->get(8)->countChildren(true) == 1 );
  }

  void testLocalNavigation() {
    Node* n1 = new Node(7.5);
    n1->make_local();
    Node* n2 = new Node(6.5);
    n2->make_nonlocal(1.0);

    Node *n3 = forest->nodes()->at(4), 
         *root = forest->local_root(),
         *n4 = forest->nodes()->at(5);

    root->remove_child(n3);
    forest->addNodeToTree(n1, root, n3, NULL);
    forest->addNodeToTree(n2, n1, NULL, NULL);

    CPPUNIT_ASSERT(n1->countChildren() == 2);
    CPPUNIT_ASSERT(n1->countChildren(true) == 1);
    CPPUNIT_ASSERT(n1->getLocalParent() == root);
    CPPUNIT_ASSERT(n3->getLocalParent() == root);

    CPPUNIT_ASSERT( (root->getLocalChild1() == n3 && root->getLocalChild2() == n4) || 
                    (root->getLocalChild1() == n4 && root->getLocalChild2() == n3) );

    CPPUNIT_ASSERT(forest->nodes()->at(0)->getLocalChild1() == NULL);
    CPPUNIT_ASSERT(forest->nodes()->at(1)->getLocalChild1() == NULL);
    CPPUNIT_ASSERT(n4->getLocalParent() == root);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestNode );

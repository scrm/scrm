#include <iostream>
#include "forest.h"
#include "model.h"
#include "node.h"
#include "random/random_generator.h"
#include "random/mersenne_twister.h"
#include "random/constant_generator.h"
#include "random/fake_generator.h"

using namespace std;

#ifndef UNITTEST
int main () {
    Model model = Model(5);

    //Different random generators.
    //MersenneTwister *rg = new MersenneTwister(5);
    //ConstantGenerator *rg = new ConstantGenerator();
    FakeRandomGenerator *rg = new FakeRandomGenerator(5);

    //Forest forest = Forest(model, rg);
    //forest.buildInitialTree();
    //double height_above = 0;
    //Node* node = new Node();
    //forest.samplePoint(true, &node, &height_above);
    //cout << node << " " << height_above << endl;
 
    cout << rg->sample() << endl;

    Forest forest = Forest(model, rg);
    forest.buildInitialTree();

    return 0;
}
#endif

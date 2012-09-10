#include <iostream>
#include "forest.h"
#include "model.h"
#include "node.h"
#include "random.h"
#include "fakerandom.h"

using namespace std;

#ifndef UNITTEST
int main () {
    Model model = Model(5);
    FakeRandomGenerator rg = FakeRandomGenerator(10);
    rg.initialize();
    cout << "rg created" << endl;
    cout << rg.sample() << endl;
    cout << rg.sampleExpo(5) << endl;
    //Forest forest = Forest(model, rg);
    //forest.buildInitialTree();
    return 0;
}
#endif

#include <iostream>
#include "forest.h"
#include "model.h"
#include "node.h"
#include "random.h"

using namespace std;

#ifndef UNITTEST
int main () {
    Model model = Model(5);
    RandomGenerator rg = RandomGenerator();
    Forest forest = Forest(model, rg);
    forest.buildInitialTree();
    return 0;
}
#endif

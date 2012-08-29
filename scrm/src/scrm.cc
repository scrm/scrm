#include <iostream>
#include "forest.h"
#include "model.h"
#include "node.h"
using namespace std;

#ifndef UNITTEST
int main () {
    Model model = Model(5);
    Forest forest = Forest(model);
    forest.buildInitialTree();
    return 0;
}
#endif

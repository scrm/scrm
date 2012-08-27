#include <iostream>
#include "node.h"
using namespace std;

#ifndef UNITTEST
int main () {
    Node node1, node2, node3;
    node1.set_height(1);
    cout << "Height" << node1.get_height() << endl;
    return 0;
}
#endif

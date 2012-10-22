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
    try {
      Model model = Model(5);

      //Different random generators.
      MersenneTwister *rg = new MersenneTwister();
      //ConstantGenerator *rg = new ConstantGenerator();
      //FakeRandomGenerator *rg = new FakeRandomGenerator(5);

      Forest forest = Forest(model, rg);
      forest.buildInitialTree();
      forest.checkTree();
      //forest.sampleNextGenealogy();

      return 0;
    }
    catch (const exception &e)
    {
      cerr << "Error: " << e.what() << endl;
      return EXIT_FAILURE;
    }
}
#endif

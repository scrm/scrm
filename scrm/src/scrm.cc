#include <iostream>
#include "forest.h"
#include "model.h"
#include "node.h"
#include "random/random_generator.h"
#include "random/mersenne_twister.h"
#include "random/constant_generator.h"
#include "random/fake_generator.h"

#ifndef UNITTEST
int main () {
    try {
      Model model = Model(10);

      //Different random generators.
      MersenneTwister *rg = new MersenneTwister();
      //ConstantGenerator *rg = new ConstantGenerator();
      //FakeRandomGenerator *rg = new FakeRandomGenerator(5);

      Forest forest = Forest(model, rg);
      forest.buildInitialTree();
      //forest.createExampleTree();
      //forest.printNodes();
      //forest.printTree();
      //forest.createExampleTree();
      
      dout << "----- NEXT -----" << std::endl;
      forest.sampleNextGenealogy();

      return 0;
    }
    catch (const exception &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}
#endif

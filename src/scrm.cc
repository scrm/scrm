#include <iostream>
#include "forest.h"
#include "random/random_generator.h"
#include "random/mersenne_twister.h"
#include "random/constant_generator.h"
#include "random/fake_generator.h"
#include <ctime>

#ifndef UNITTEST
int main () {
    try {
      time_t start_time = time(0);
      Model model = Model(5);

      //Different random generators.
      MersenneTwister *rg = new MersenneTwister();
      //ConstantGenerator *rg = new ConstantGenerator();
      //FakeRandomGenerator *rg = new FakeRandomGenerator();

      Forest forest = Forest(model, rg);
      forest.buildInitialTree();
       
      for (size_t i=1; i <= 25; ++i) {
        forest.set_current_base(i);    
        forest.sampleNextGenealogy();
      }

      std::cout << forest.getNodes()->size() << std::endl;

      time_t end_time = time(0);
      std::cout << "Simulation took about " << end_time - start_time 
                << " second(s)" << std::endl;

      delete rg;
      return 0;
    }
    catch (const exception &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}
#endif

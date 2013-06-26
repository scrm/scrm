#ifndef scrm_param
#define scrm_param

#include <iostream>
#include <iomanip>      
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>  
#include <stdio.h>
#include <stdexcept>
#include <boost/lexical_cast.hpp> 

//#include "mtrand.h"
#include "model.h"

class Param {
 public:
  bool seg_bool;
  int total_mut;
  int random_seed;
  std::string treefile;
  
  bool tmrca_bool;
  std::string tmrca_NAME;
  
  void init();

  Param() : argc_(0), argv_(NULL) { };
  Param(int argc, char *argv[]) : argc_(argc), argv_(argv) { }

  Model parse();
  
  void print_param();
  
  bool log_bool;
  std::string log_NAME;
  
  void log_param();
   
 private:
  const int argc_;
  char * const* argv_;
};

void print_help();
void print_example();
void print_option();	

template<class T>
T readInput(char input[])
{
  return boost::lexical_cast<T>(input);
}

#endif

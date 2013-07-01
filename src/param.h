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
#include "model.h"

class Param {
 public:
  friend std::ostream& operator<< (std::ostream& stream, const Param& param) {
    stream << "scrm ";
    for (size_t i = 1; i < param.argc_; ++i) {
      stream << param.argv_[i] << " ";
    }
    return stream;
  };

  int random_seed;
  std::string treefile;
  
  std::string tmrca_NAME;
  
  void init();

  bool seg_bool() const { return seg_bool_; }
  bool tmrca_bool() const { return tmrca_bool_; }

  void set_seg_bool(const bool &seg_bool) { seg_bool_ = seg_bool; } 
  void set_tmrca_bool(const bool &tmrca_bool) { tmrca_bool_ = tmrca_bool; } 

  Param() : argc_(0), argv_(NULL) { };
  Param(int argc, char *argv[]) : argc_(argc), argv_(argv) { }

  Model* parse();
  
  void print_param();
  
  bool tree_bool;
  bool log_bool;
  std::string log_NAME;
  
  void log_param();
   
 private:
  const int argc_;
  char * const* argv_;

  bool tmrca_bool_;
  bool seg_bool_; 
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

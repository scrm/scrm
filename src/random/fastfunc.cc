#include "fastfunc.h"
#include "mersenne_twister.h"

float* FastFunc::build_fastlog_table() {
  float * table = (float*) malloc( 1025*sizeof(float) );
  for (int index=0; index<1025; index++) {
    double dlog1 = log( 1.0 + (index) / 1024. );
    double dlog2 = log( 1.0 + (index+1.0) / 1024. );
    double middiff1 = log( 1.0 + (index+0.5)/1024. ) - (dlog1 + (dlog2-dlog1)*0.5);
    table[index] = (float)(dlog1 + middiff1/2);
  }
  return table;
}

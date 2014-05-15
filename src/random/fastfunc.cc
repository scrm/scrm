/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013, 2014 Paul R. Staab, Sha (Joe) Zhu and Gerton Lunter
 * 
 * This file is part of scrm.
 * 
 * scrm is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "fastfunc.h"

float* FastFunc::build_fastlog_table() {
  float * table = (float*) malloc( 1025*sizeof(float) );
  double prevx = 1.0;
  double prevy = 0.0;
  for (int index=0; index<1025; index++) {
    /* 
     * Previous algorithm.  This has good maximum deviation from the
     * true logarithm, but results in log(1.0) being just slightly >0,
     * causing a problem in the sampling algorithms.
     *
    double dlog1 = log( 1.0 + (index) / 1024. );
    double dlog2 = log( 1.0 + (index+1.0) / 1024. );
    double middiff1 = log( 1.0 + (index+0.5)/1024. ) - (dlog1 + (dlog2-dlog1)*0.5);
    table[index] = (float)(dlog1 + middiff1/2);
     *
     */

    // calculate x coordinate at which linear approximation is exactly equal
    // to the true logarithm
    double curx = 1.0 + (index+5.0/6.0) / 1024;
    if (index == 1023)
      curx = 1.0 + (index+1.0) / 1024;

    // calculate true logarithm at the point
    double cury = log( curx );

    // calculate the linear approximation at the next join point
    double targetx = 1.0 + (index+1.0)/1024;
    double targety = prevy + (targetx-prevx)*(cury-prevy)/(curx-prevx);

    // store previous linear approximation in table, and update prevx/y
    table[index] = (float)(prevy);
    prevx = targetx;
    prevy = targety;
  }
  return table;
}

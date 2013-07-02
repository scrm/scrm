/*
 * scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.
 * 
 * Copyright (C) 2013 Paul R. Staab, Sha (Joe) Zhu and Gerton Lunter
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

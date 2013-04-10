/*
    ICSI benchmark tool to test the speed up on your system
    Test using a set of random numbers of different logarithm implementations.

    Version 0.6 beta
    Build date: June 17st, 2007

    Copyright (C) 2007 International Computer Science Institute
    1947 Center Street. Suite 600
    Berkeley, CA 94704
    
    Contact information: 
         Oriol Vinyals	vinyals@icsi.berkeley.edu
         Gerald Friedland 	fractor@icsi.berkeley.edu
 
    Acknowledgements:
    Thanks to Harrison Ainsworth (hxa7241@gmail.com) for his idea that
    doubled the accuracy (icsilog v2).

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "icsilog.h"


double currenttime();
time_t time (time_t * timer);
inline float fast_log (float val);
 

int main(int argc, char *argv[])
{
   int iters,N,i;
   float randomrange;
   float *result1,*result2,*result3,*result4,*result5,*randomv,*lograndom,*LOOKUP_TABLE, *LOOKUP_TABLE2;
   double error1,error2,error3,error4,error5,start,end,time1,time2,time3,time4,time5;

   if(argc!=4)
   {
	printf("Usage: %s NumIters BitsTable MaxRange\n",argv[0]);
	printf("\nNumIters: Total number of random numbers to be tested.\n");
	printf("BitsTable: Number of bits of the mantissa. Table size will have 2^BitsTable positions.\n");
	printf("MaxRange: Maximum value of the generated random numbers.\n");
        exit(1);
   } 
   srand(time(0)); /*change seed*/

   /*INITIALIZE variables*/
   iters = atoi(argv[1]);
   N = atoi(argv[2]);
   if(N>23 || N < 0)
   {
    	printf("Input error: Second parameter (BitsTable) must be between 0 and 23\n");
	exit(1);
   }
   randomrange = atof(argv[3]);
   result1 = (float*) malloc(iters*sizeof(float));
   result2 = (float*) malloc(iters*sizeof(float));
   result3 = (float*) malloc(iters*sizeof(float));
   result4 = (float*) malloc(iters*sizeof(float));
   result5 = (float*) malloc(iters*sizeof(float));
   randomv = (float*) malloc(iters*sizeof(float));
   lograndom = (float*) malloc(iters*sizeof(float));
   error1 = 0;
   error2 = 0;
   error3 = 0;
   error4 = 0;
   error5 = 0;
 
   /*Construct the lookup table*/
   LOOKUP_TABLE = (float*) malloc(((int) pow(2,N))*sizeof(float));
   LOOKUP_TABLE2 = (float*) malloc(((int) pow(2,N))*sizeof(float));
   fill_icsi_log_table(N,LOOKUP_TABLE);
   fill_icsi_log_table2(N,LOOKUP_TABLE2);

   /*Save the exact (double precision) value of the log*/
   for(i=0;i<iters;i++)
   {
	randomv[i] = (float)rand() / RAND_MAX * randomrange;
        lograndom[i] = log(randomv[i]);
   }

   /*Experiment starts*/
   start = currenttime();
   for(i=0;i<iters;i++)
   {
   	result1[i] = fast_log(randomv[i]);
   }
   end = currenttime();
   time1 = end - start;
   start = currenttime();
   for(i=0;i<iters;i++)
   {
   	result2[i] = icsi_log(randomv[i],LOOKUP_TABLE,N);
   }
   end = currenttime();
   time2 = end - start;
   start = currenttime();
   for(i=0;i<iters;i++)
   {
   	result3[i] = logf(randomv[i]);
   }
   end = currenttime();
   time3 = end - start;

   start = currenttime();
   for(i=0;i<iters;i++)
   {
   	result4[i] = log(randomv[i]);
   }
   end = currenttime();

   time4 = end - start;

   start = currenttime();
   for(i=0;i<iters;i++)
   {
   	result5[i] = icsi_log_v2(randomv[i],LOOKUP_TABLE2,N);
   }
   end = currenttime();

   time5 = end - start;

   /*Compute the errors*/
   for(i=0;i<iters;i++)
   {
   	error1 += pow(result1[i]-lograndom[i],2.0);
   	error2 += pow(result2[i]-lograndom[i],2.0);
   	error3 += pow(result3[i]-lograndom[i],2.0);
   	error4 += pow(result4[i]-lograndom[i],2.0);
   	error5 += pow(result5[i]-lograndom[i],2.0);
   }
   error1=sqrt(1.0/iters*(error1));
   error2=sqrt(1.0/iters*(error2));
   error3=sqrt(1.0/iters*(error3));
   error4=sqrt(1.0/iters*(error4));
   error5=sqrt(1.0/iters*(error5));

   /*Report results*/
   printf("#### Test of several single precision logarithm functions ####\n");
   printf("----------------------------------------------------------------------------------------------------------\n");
   printf("Experiment with %d samples and LUT of %d positions\n",iters,(int) pow(2,N));
   printf("----------------------------------------------------------------------------------------------------------\n");
   printf("log(): Time: %f secs. Total error: %e. SpeedUp with respect to log(): %f\n",time4,error4,time4/time4);   
   printf("----------------------------------------------------------------------------------------------------------\n");
   printf("logf(): Time: %f secs. Total error: %e. SpeedUp with respect to log(): %f\n",time3,error3,time4/time3); 
   printf("----------------------------------------------------------------------------------------------------------\n");
   printf("Fast Log: Time: %f secs. Total error: %e. SpeedUp with respect to log(): %f\n",time1,error1,time4/time1);  
   printf("----------------------------------------------------------------------------------------------------------\n");
   printf("ICSI Log: Time: %f secs. Total error: %e. SpeedUp with respect to log(): %f\n",time2,error2,time4/time2);
   printf("----------------------------------------------------------------------------------------------------------\n");
   printf("ICSI Log v2: Time: %f secs. Total error: %e. SpeedUp with respect to log(): %f\n",time5,error5,time4/time5);
   printf("----------------------------------------------------------------------------------------------------------\n");
   return(0);
}

/*Taylor (deg 3) implementation of the log: http://www.flipcode.com/cgi-bin/fcarticles.cgi?show=63828*/
inline float fast_log (float val)
{
   register int *const     exp_ptr = ((int*)&val);
   register int            x = *exp_ptr;
   register const int      log_2 = ((x >> 23) & 255) - 128;
   x &= ~(255 << 23);
   x += 127 << 23;
   *exp_ptr = x;

   val = ((-1.0f/3) * val + 2) * val - 2.0f/3;

   return ((val + log_2)* 0.69314718);
}

/*timing function*/
double currenttime()
{
  struct timeval curtime;
  gettimeofday(&curtime,NULL); 
  return curtime.tv_sec+(curtime.tv_usec/1000000.0);
}

/*  stats.c : Calculates mean, standard deviation for a list of numbers.
   Quantiles are also calculated if quantiles are specified by command 
line arguments. For example,  stats 0.05  0.5  0.95 <datafile
would output the mean, standard deviation (estimated from sample) and 
estimates of the  0.5, 0.5 and 0.95th quantile.  
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

main( int argc, char *argv[])
{

	double *vec, x, s, *percentiles ;
	int  c, vecl = 1000 ;
	int i, index ;
	double p;
	int order(int, double *) ;

	percentiles = (double *)malloc( (unsigned)argc*sizeof(double) );
	vec = (double *)malloc( (unsigned)vecl*sizeof( double) ) ;

	for( i=1; i<argc; i++){
	  percentiles[i] = atof( argv[i] ) ;
	  }
	c = 0 ;
	x = s= 0.0 ;

	while( scanf(" %lf",vec+c) != EOF ) {
	   x += vec[c];
	  s += vec[c]*vec[c] ;
           c++;
	 if( c >= vecl ) {
	    vecl += 1000;
	    vec = (double *)realloc( vec, (unsigned)vecl*sizeof(double) ) ;
	   }
	}	
	 order( c, vec);
	x /= c ;
	s /= c ;
	s -= x*x ;
	s = sqrt( s*c/(c-1.0) ) ;
	printf("%lf\tsd:\t%lf\tn:\t%d", x, s, c);

	for( i=1; i<argc; i++){
	  index = percentiles[i]*c + 0.5  ;
	  p = vec[index-1]*(index + 0.5 - percentiles[i]*c) 
		 + vec[index]*(percentiles[i]*c+0.5 -index) ;
	   printf("\t%5.3lf", percentiles[i]);
	  if( index < 1 ) printf("\t-");
	   else printf("\t%lf", p);
	}
	printf("\n");
}	


       int
order(int n, double *pbuf)
{
        int gap, i, j;
        double temp;

        for( gap= n/2; gap>0; gap /= 2)
           for( i=gap; i<n; i++)
                for( j=i-gap; j>=0 && pbuf[j]>pbuf[j+gap]; j -=gap) {
                   temp = pbuf[j];
                   pbuf[j] = pbuf[j+gap];
                   pbuf[j+gap] = temp;
                   }
}



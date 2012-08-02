/* *********************** Probability of s segregating sites (no recombination) ********* 

	The subroutine dist_ss() calculates the probabilities of number of segregating sites, using
a simple recursion (Hudson, RR in Oxford Surveys in Evol. Biol. V.7).  


********************************************************************************** */

#include <stdio.h>
  #include <stdlib.h>

main(int argc, char *argv[])
{
	double theta, *p1, cum_prob, ess, varss, e_ss(), var_ss(), sqrt(), r,
		 varr, varssr() ;
	int nsam, mmax, m, i  ;

	if( argc < 3) {
	    printf("Usage: dists3 nsam s theta1 theta2 ... \n");
		exit(0);
		}
	nsam = atoi( argv[1] ) ;
	mmax = atoi( argv[2] );

		p1 = (double *)malloc( (unsigned)((mmax+1)*sizeof(double)) ) ;

	for( i= 3; i<  argc; i++){
	   dist_ss(theta =atof(argv[i]),nsam,mmax,p1);
	   printf("%lf\t%lf\n",theta, p1[mmax] );
	   }
	
}
	
	
	int
dist_ss(theta,nsam,m,p1)
	double theta, *p1;
	int nsam, m;
{
	double  *p2, qjm, *ptmp, c, sum, *pin, *pnew ;
	int n, i, j,  k ;

	pin = p1 ;
	p2 = (double *)malloc( (unsigned)((m+1)*sizeof(double)) ) ;
	pnew = p2;
	
	c = theta/(theta + 1.) ;
	p2[0] =1./(theta + 1.) ;
	for(i=1;i<=m; i++)
		p2[i] = c*p2[i-1];
	
	if( nsam > 2 ) {
	 for( n=3; n<=nsam; n++){
	  for( k=0; k<=m ; k++) {
	  	p1[k] = 0.0;	
	 	qjm = (n-1.)/(theta+n-1.) ;
	    for( j=0; j<=k ; j++) {
	   	  p1[k] += qjm*p2[k-j] ;
	   	  qjm *= theta/(theta+n-1.) ;
	   	  }
	    }
	  ptmp = p2 ;
	  p2 = p1;
	  p1 = ptmp;
	  }
	 }
	 for( i=0; i<=m; i++) pin[i] = p2[i] ;
	 free(pnew);
}

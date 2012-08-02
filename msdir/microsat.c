/*  microsat.c:   This function takes ms output and converts it is microsatellite
  length variation data.  The output has on each line the set of lengths
  of the nsam individuals (relative to the ancestral length).
  Example usage:   ms 10 5 -t 4.0 | microsat > msat.dat
  To compile:  gcc -o microsat microsat.c rand1.c -lm   
*/

#include <stdio.h>
#include <stdlib.h>


int maxsites = 1000 ;
double ran1() ;

main(argc,argv)
	int argc;
	char *argv[];
{
	int nsam, j ,nsites, i,  howmany  ;
	char **list, **cmatrix(), allele,na, line[1001]  ;
	FILE *pf, *fopen(), *pfin ;
	double *posit   ;
	int   segsites, count  , nadv, probflag  ;
	double prob ;
	char dum[20], astr[100] ;
	int *nrepeats, step, ind ;

/* read in first two lines of output  (parameters and seed) */
  pfin = stdin ;
  fgets( line, 1000, pfin);
  sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
  fgets( line, 1000, pfin);

	if( argc > 1 ) { 
	   nadv = atoi( argv[1] ) ; 
	}

  list = cmatrix(nsam,maxsites+1);
  posit = (double *)malloc( maxsites*sizeof( double ) ) ;
  nrepeats = (int *)malloc(nsam*sizeof(int) );

  count=0;
	probflag = 0 ;
while( howmany-count++ ) {

/* read in a sample */
  do {
     if( fgets( line, 1000, pfin) == NULL ) exit(0);
  }while ( line[0] != '/' );
 
  fscanf(pfin,"  segsites: %d", &segsites );
  if( segsites >= maxsites){
	maxsites = segsites + 10 ;
	posit = (double *)realloc( posit, maxsites*sizeof( double) ) ;
        biggerlist(nsam,maxsites, list) ;
        }
   if( segsites > 0) {
	fscanf(pfin," %s", astr);
	if( astr[1] == 'r' ){
	   fscanf(pfin," %lf", &prob ) ;
	   probflag = 1;
	   fscanf(pfin," %*s");
	}
	for( i=0; i<segsites ; i++) fscanf(pfin," %lf",posit+i) ;
	for( i=0; i<nsam;i++) fscanf(pfin," %s", list[i] );
	}
/* analyse sample ( do stuff with segsites and list) */
   for( ind = 0; ind < nsam; ind++) nrepeats[ind] = 0 ;
   for( i = 0; i< segsites; i++){
     if( ran1() < .5) step = -1 ;
     else step = 1 ;
     for( ind = 0; ind < nsam; ind++)
       if( list[ind][i] == '1' ) nrepeats[ind] += step ;
   }
   for( ind=0; ind < nsam-1; ind++) printf("%d\t",nrepeats[ind] );
   printf("%d",nrepeats[nsam-1]);
   printf("\t%s",line+2 );
  }
}

	

/* allocates space for gametes (character strings) */
	char **
cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned)( nsam*sizeof( char* )) ) ) )
	   perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) (len*sizeof( char )) )))
			perror("alloc error in cmatric. 2");
		}
	return( m );
}

        int
biggerlist(nsam, nmax, list )
        int nsam ;
        unsigned nmax ;
        char ** list ;
{
        int i;

        maxsites = nmax  ;
        for( i=0; i<nsam; i++){
           list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
           if( list[i] == NULL ) perror( "realloc error. bigger");
           }
}                        


# In this file are scripts that use a simple rejection algorithm
# to obtain draws from the posterior distribution of a demographic
# parameter. 

# The problem:  With observed values of Tajima's D from 40 unlinked loci,
# make inferences about the time of expansion of a population. 

# Solution:  Summarize the data by the average Tajima's D (over the 40 loci) and
# use an approximate rejection method to obtain a posterior distribution of the
# time of expansion given the observed average value of Tajima's D.

# Assumptions:  A single panmictic population with theta = 30 for all loci, no
# recombination within loci, sample size is 50 at all loci, and the expansion was an
# instantaneous increase by a factor of 20, that is, the ancient population size 
# divided by the current size is 0.05.  The prior distribution of the time of expansion
# is assumed uniformly distributed in the interval (0., tmax).

# The following R scripts generate independent draws from this posterior distribution.
# The function mtajd1 calls ms and sstats which must be
# in the working directory. ( sstats is just the name I give to the executable based on
# sample_stats.c  It can be generated as follows:  gcc -o sstats sample_stats.c tajd.c -lm  )
# The function posttimes generates the
# independent draws.  An example of the use of the function posttimes
# appears near the bottow of this file.

# The following function, mtajd1 generates independent polymorphism data sets,
# that are interpreted as independent loci, each with theta=30, and rho=0.0. The two
# arguments to the function are texpansion, the time of expansion, and nloci, the number of
# loci. The data are generated with calls to ms, the output of ms is piped to sstats, 
# which calculates Tajima's D and other stats for each locus. The Tajima's D values are cut
# out of the output with the "cut" function and assigned to the local R variable tajd.  
# The mean value is then calculated and returned. This routine doesn't use the "tbs" 
# command line argument. 


mtajd1 <- function( texpansion, nloci=20, nsam=50, theta=30.) {
	mscall <- paste( "./ms ",nsam,nloci," -t ",theta, "  -eN ",texpansion, " .05 | ./sstats | cut -f 6 ")
	tajd <- system( mscall, intern=TRUE)
  mtajd <- mean( as.numeric( tajd) )
  mtajd
  }

# The following function returns ndraws from the posterior distribution of the 
# time of expansion conditional on the mean TD being within eps of an 
# observed value.

posttimes <- function( ndraws = 100, mtd_obs, tmax=0.2,eps=.1, nloci=20)  {
 postv <- numeric()
 mycount <- 0
 for( i in 1:ndraws ){
 repeat {
 	rtime <- runif(1,0,tmax)   # generate a random time of expansion
 	mycount <- mycount +1
    td <- mtajd1(rtime,nloci)   # generate a random mean tajima's D given texp
    if(  abs( mtd_obs -td) < eps  ) {    # if within epsilon, accept
    	postv[i] <- rtime
    	break
    	}
    }
 }
 print( ndraws/mycount )   # prints out the fraction of iterations accepted
 postv
}

# An example application of posttimes, showing the posterior distribution
# of texpansion given observed mean Taj''s D is -1.5 (for 40 loci).  This takes
# few minutes.
 
  mypost <- posttimes(ndraws=400,mtd_obs=-1.5, nloci=40)
  plot(density(mypost), xlab="time of expansion" )
 rug(mypost)


# The following bit of code plots estimated mean Taj''s D  as a function of texpansion.  This 
# is here to help one understand the results from the above example. (The plot produced 
# shows that mean Tajima's D is around -1.5 for texpansion = 0.009 and for texpansion = 0.15 .)
# This takes a minute or two.

 texp <- c( (1:15)/1000, (2:20)/100 )
 mymeans <- sapply(texp,mtajd1,1000)
 plot(texp,mymeans, xlab="time of expansion", ylab="mean Tajima's D")


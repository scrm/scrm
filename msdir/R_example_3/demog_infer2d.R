# A simple example of a rejection algorithm for demographic inference of two parameters.
# This function returns ndraws from the posterior distribution of the 
# time of expansion and f (the ratio of  ancestral population size to present population
# size) conditional on the mean TD being within eps of an 
# observed value. (In this example, the mean is over 20 unlinked loci.
# The prior distribution of the time is uniform (0,tmax)
# and the prior for f is uniform (0,fmax).  The population size change is instantaneous.
# In this example, we assume that the theta's for the set of loci are uniformly 
# distributed on the interval (5,40), and the rho's are uniformly distributed on (0,10).
# The executables ms and sstats must be in the working directory.
# sstats is the executable based on sample_stats.c, and is produced with
# gcc -o sstats sample_stats.c tajd.c -lm  

# Here is the function that generates the draws from the posterior distribution.
# It calls mtajd() which follows.

posttimes <- function( mtd_obs, ndraws = 100, tmax=0.2,eps=.1, nloci=20, fmax=0.2)  {
 posttexp <- numeric()
 postf <- numeric()
 mycount <- 0
 for( i in 1:ndraws ){
 repeat {
 	rtime <- runif(1,0,tmax)  # a draw from the prior distribution of texpansion
 	rf <- runif(1,max=fmax)   # from the prior for f
 	mycount <- mycount +1
    td <- mtajd(rtime,rf,nloci)  # returns a random realization of mean tD
    if(  abs( mtd_obs -td) < eps  ) {  # If close enough to observed, accept
    	posttexp[i] <- rtime
    	postf[i] <- rf
    	break
    	}
    }
 }
 print( ndraws/mycount )
 cbind(posttexp,postf)
}

# this function calls ms and sstats to obtain a mean Tajima s D 
#  value of a set of unlinked loci, under a model with recent population
# expansion.  The time of expansion is in units of 4No generations.  
# In this example, each loci has a random theta and rho values. nsam=50 .
# Theta and rho are random. rf is the ratio of ancient to present population size.

mtajd <- function( texpansion,rf, nloci=20, nsam=50) {
        # theta and rho values are randomly assigned from uniform distrib.
	thetavalues <- runif(nloci,min=5.0, max=40.)  # theta uniform on (5,40) 
	rhovalues <- runif(nloci,max=10.)             # rho uniform on (0,10)
	paramvalues <- cbind(thetavalues,rhovalues)
	mscall <- paste( "| ./ms ", nsam ,nloci," -t tbs -r tbs 1000 -eN ",texpansion,rf, " | ./sstats | cut -f 6 >tajd.out")
	write(t(paramvalues),file=mscall)
   tajd <- scan("tajd.out", quiet=TRUE)
   mtajd <- mean(  tajd )
   mtajd
  } 

#  Example using posttimes, and plotting results. This takes ten minutes or so
#   on my little old laptop.

mypost <- posttimes(-1.5,ndraws=2000, nloci=20)
 # to plot an estimate of 2-dim density:
  library(KernSmooth)
  est <- bkde2D(mypost,bandwidth=c(.01,.01))
  contour(est$x1,est$x2,est$fhat, xlab="t-expansion",ylab="f") 
   # or for a plot of the points, uncomment the following:
# plot(mypost, xlab="texpansion",ylab="f")


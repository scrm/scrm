# An example of the Metropolis-Hastings MCMC without likelihoods
# method ( Marjoram, P., Molitor, J., Plagnol, V. & Tavar´e, S. Markov
# chain Monte Carlo without likelihoods. Proc. Nat. Acad. Sci. USA 100, 15324–15328 (2003). )

# This function returns ndraws from an mcmc process to estimate the posterior
# distribution of the
# time of expansion and f (the ratio of the ancestral to present 
# population size) conditional on the mean TD over nloci 
# being within eps of an 
# observed value.  The prior distribution of the time is uniform (0,tmax)
# and the prior for f is uniform (0,fmax).  The proposal distribution is
# symmetric, and the priors are uniform,
#  thus the metropolis-hastings step is simple.
# In this simple version there is no burn-in and the sampling interval is one.
# The executables ms, sstats must be in the working directory.
# We assume the population size change is instantaneous, we have mean Tajima's 
# D from 40 loci, thetas for the loci are from a uniform on (5,40), and the rho's
# are uniform from (0,10).

demog_mcmc <- function( mtd_obs, niters = 100, tmax=0.2,eps=.1, nloci=20, fmax=0.2)  {
 posttexp <- numeric()
 postf <- numeric()
 texpold <- tmax/2.0  # starting value of texp
 fold <- fmax/2.0     # starting value of f
 mycount <- 0
 for( i in 1:niters ){
 	texpnew <- abs( texpold + rnorm(1,sd=tmax/15.) )  # sd is important for how this process mixes
 	if( texpnew > tmax ) texpnew <- 2*tmax - texpnew   # reflection at boundary
 						# The above line can fail if texpnew is bigger than 2*tmax 
 	fnew <- abs( fold + rnorm(1,sd=fmax/15.) )
 	if( fnew > fmax) fnew <- 2*fmax - fnew             # reflection at boundary
    td <- mtajd(texpnew,fnew,nloci) # returns a random realization, for mean tD
    if(  abs( mtd_obs -td) < eps  ) {  # If close enough to observed, accept
    	texpold <- texpnew
    	fold <- fnew
 	    mycount <- mycount +1
    }
    posttexp[i] <- texpold
    postf[i] <- fold
 }
 print( mycount/niters )
 cbind(posttexp,postf)
}

# this function calls ms and sstats to obtain a mean Tajima s D 
#  value of a set of unlinked loci, under a model with recent population
# expansion.  The time of expansion is in units of 4No generations.  
# In this example, each loci has a random theta and rho values. nsam=50 .

mtajd <- function( texpansion,rf, nloci=20) {
        # theta and rho values are randomly assigned from uniform distrib.
	thetavalues <- runif(nloci,min=5.0, max=40.) 
	rhovalues <- runif(nloci,max=10.)
	paramvalues <- cbind(thetavalues,rhovalues)
	# With the following two lines, the theta and rho values are piped to mstbs 
	mscall <- paste( "| ./ms 50 ",nloci," -t tbs -r tbs 1000 -eN ",texpansion,rf, " | ./sstats | cut -f 6 >tajd.out") 
	write(t(paramvalues),file=mscall)
   tajd <- scan("tajd.out", quiet=TRUE)
   mtajd <- mean(  tajd )
   mtajd
  }

# To make a plot of mean taj s D as a function of texpansion, for f=0.05
 
 texp <- c( (1:10)/1000, (1:20)/100 )
 mymeans <- sapply(texp,mtajd,0.05,1000)
 plot(texp,mymeans,xlab="texp", ylab="T_D")

#  Example using demog_mcmc, and plotting results

mymcmc <- demog_mcmc(-1.5,niters=4000,nloci=40)  # With 4000 iterations this takes a while (maybe 10 minutes) but makes a nice plot.
plot(mymcmc, xlab="texpansion",ylab="f" )   # plot the points, or alternatively, plot contour lines as follows:
library(KernSmooth)
est <- bkde2D(mymcmc,bandwidth=c(.01,.01))
contour(est$x1,est$x2,est$fhat, xlab="texpansion",ylab="f" )  


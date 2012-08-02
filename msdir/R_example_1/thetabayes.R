# Inference about theta( 4Nu) from an observation of the number of segregating 
#     sites at a locus with recombination.
#
#  The following script will generate independent draws from the
# posterior distribution of theta conditional on the observed 
# number of segregating sites is equal to 14 (nsam=10) and assuming the prior
# distribution of theta is uniform on the interval (0,20) and the recombination
# parameter (rho=4Nr) equals 20. The method is a simple rejection scheme.  
# From these draws from the posterior distribution one can estimate the density
# of the posterior distribution and plot this distribution.  The script plots this
# estimated distribution.  The script below also generates draws from the posterior
# distribution under the
# assumption that the recombination rate is zero.  For this case the probabiility
# of 14 segregating sites can be calculated and compared to the above densities.
# This is also done with the following script.
#
# The executables ms and dist3 need to be in the workng directory.
#  These are compiled as follows:
#  ./clms
#  gcc -o dist3 dist3.c -lm
#
#  The pdf file, thetabayesfig.pdf shows the plot produced by this script.

# These take a minute or two.
#  For rho=20.  :
myunif <- runif(30000,max=20.0)   # draws from the prior distribution of theta
 write(myunif,file="|./ms 10 30000 -t tbs -r 20. 1000 | grep 'segsites: 14$' -B 1 | grep '//' | cut -f 2  >msout")
 thetas <- scan("msout")    # posterior theta''s 
 plot(density(thetas), xlim=c(0,20), xlab="theta" )
 rug(thetas)

#  For rho=0.  : 
myunif <- runif(30000,max=20.0)
 write(myunif,file="|./ms 10 30000 -t tbs -r 0. 1000 | grep 'segsites: 14$' -B 1 | grep '//' | cut -f 2  >msout")
 thetas <- scan("msout")
 lines(density(thetas), col="green" )

# The following calculates P(S=14| theta) for theta = 1, 2, 3, ..., 20 under
#  the assumption that rho=0, and plots
# these values (normalized by the sum of all of them.)  Uses the program dist3.c 
# Plots the results on the above plot.

mycall2 <- paste(c("./dist3 10 14", 1:20, " >mylik.out"), collapse=" " )
system(mycall2)
mylik1 <- scan("mylik.out")
dim(mylik1) <- c(2,20)
mylik2 <- t( mylik1)
points( mylik2[,1],mylik2[,2]/sum(mylik2[,2]) )

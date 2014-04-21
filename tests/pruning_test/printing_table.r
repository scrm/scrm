rm(list=ls())
library(matrixStats)

ms_     = read.table( "ms_moment" ,         sep = ",") 
scrm_   = read.table( "scrm_moment" ,       sep = ",") 
scrm_p_ = read.table( "scrmprune_moment" , sep = ",") 

meanMat = matrix(0,4,3)
meanMat[,1]=colMeans(ms_)
meanMat[,2]=colMeans(scrm_)
meanMat[,3]=colMeans(scrm_p_)

sdMat = matrix(0,4,3)
sdMat[,1]=colSds(ms_)
sdMat[,2]=colSds(scrm_)
sdMat[,3]=colSds(scrm_p_)


printmat = matrix(0,4,6)
printmat[,1]=colMeans(ms_)
printmat[,2]=colSds(ms_)
printmat[,3]=colMeans(scrm_)
printmat[,4]=colSds(scrm_)
printmat[,5]=colMeans(scrm_p_)
printmat[,6]=colSds(scrm_p_)
print(format(printmat,digits=2))


ms     = read.table( "ms_stat"   ) 
scrm   = read.table( "scrm_stat" ) 
scrm_p = read.table( "scrmprune_stat") 

printmat2 = matrix(0, 5, 6)
printmat2[,1]=colMeans(ms[,c(2,4,6,8,10)])
printmat2[,2]=colSds(ms[,c(2,4,6,8,10)])
printmat2[,3]=colMeans(scrm[,c(2,4,6,8,10)])
printmat2[,4]=colSds(scrm[,c(2,4,6,8,10)])
printmat2[,5]=colMeans(scrm_p[,c(2,4,6,8,10)])
printmat2[,6]=colSds(scrm_p[,c(2,4,6,8,10)])

print(format(printmat2,digits=3))

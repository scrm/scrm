rm(list=ls())
library(dgof) # overwrite ks.test() function, allow discrete case. However, it does not make any difference in our case.
tmrca_moment <- function ( ms, scrm, s10000, s50000){
    ms_scrm_test = list ( ks.test( ms$V1, scrm$V1 ),
                          ks.test( ms$V2, scrm$V2 ),
                          ks.test( ms$V3, scrm$V3 ),
                          ks.test( ms$V4, scrm$V4 ) )
    
    ms_s10000_test = list ( ks.test( ms$V1, s10000$V1 ),
                           ks.test( ms$V2, s10000$V2 ),
                           ks.test( ms$V3, s10000$V3 ),
                           ks.test( ms$V4, s10000$V4 ) )    
    
    ms_s50000_test = list (ks.test( ms$V1, s50000$V1 ),
                           ks.test( ms$V2, s50000$V2 ),
                           ks.test( ms$V3, s50000$V3 ),
                           ks.test( ms$V4, s50000$V4 ) )

    titles = c("1st", "2nd", "3rd", "4th")
    
    pdf ( "TMRCA_moment_KStest.pdf" )
    par(mfrow=c(2,2))
    for ( i in c(1:4) ) {
        ms_i    = ms[,i]
        scrm_i  = scrm[,i]
        scrm_p1 = s10000[,i]
        scrm_p5 = s50000[,i]

        xrange = range ( c(ms_i, scrm_i, scrm_p1, scrm_p5) )

        plot( ecdf(ms_i), xlim = xrange, col="red", main = paste( titles[i], "Moment" ) )
        plot( ecdf(scrm_i), add=TRUE, lty="dashed", col="blue")
        plot( ecdf(scrm_p1), add=TRUE, lty="dashed", col="green")
        plot( ecdf(scrm_p5), add=TRUE, lty="dashed", col="black")

	legend( "bottomright", c( paste( "p-value = ",      format(ms_scrm_test[[i]]$p.value,   digits = 2), sep = ""),
				  paste( "p-value = ",      format(ms_s50000_test[[i]]$p.value, digits = 2), sep = ""),
                                  paste( "p-value = ",      format(ms_s10000_test[[i]]$p.value, digits = 2), sep = "") ) , 
				col=c("blue", "black", "green"), pch = 16)
        legend( "topleft" , c( "window 10000 ", "window 50000",  "scrm exact", "ms"), col=c( "green", "black", "blue", "red" ), pch=16)
    }
    dev.off()
}

ms_     = read.table( "ms_moment" ,  colClasses = "numeric",    sep = ",") 
scrm_   = read.table( "scrm_moment" , colClasses = "numeric",   sep = ",") 
#scrm_p0   = read.table( "scrmprune0_moment" , colClasses = "numeric",   sep = ",") 
scrm_p10000_ = read.table( "scrmprune10000_moment" ,colClasses = "numeric",  sep = ",") 
scrm_p50000_ = read.table( "scrmprune50000_moment" ,colClasses = "numeric",  sep = ",")

tmrca_moment ( ms_, scrm_, scrm_p10000_, scrm_p50000_)

#cat(paste(${nsam},${r},\"|\",
#format(m1test\$statistic,digits=4),format(m1test\$p.value,scientific = TRUE),\"|\",
#format(m2test\$statistic,digits=4),format(m2test\$p.value,scientific = TRUE),\"|\",
#format(m3test\$statistic,digits=4),format(m3test\$p.value,scientific = TRUE),\"|\",
#format(m4test\$statistic,digits=4),format(m4test\$p.value,scientific = TRUE),
#sep=\"\t\"),file=\"${compareTMRCA}\",append=TRUE);cat(\"\n\",file=\"${compareTMRCA}\",append=TRUE);

rm( list=ls() )

continuous <- function (ms, scrm, scrm_p10000, scrm_p50000, title_i){
    test1 = ks.test( ms, scrm )
    test2 = ks.test( ms, scrm_p50000)
    test3 = ks.test( ms, scrm_p10000)
    xrange = range ( c(ms, scrm, scrm_p10000, scrm_p50000) )
    pdf ( paste( title_i,".pdf", sep = "" ) )
    plot( ecdf(ms), xlim = xrange, col="red", main = paste( title_i ) , xlab = title_i)
    plot( ecdf(scrm), add=TRUE, lty="dashed", col="blue")
    plot( ecdf(scrm_p10000), add=TRUE, lty="dashed", col="green")
    plot( ecdf(scrm_p50000), add=TRUE, lty="dashed", col="black")
    legend( "topleft" , c( "window 10000", "window 50000", "scrm", "ms"), col=c( "green", "black",  "blue", "red" ), pch=16)
    legend( "bottomright", c( paste( "p-value = ",      format(test1$p.value,   digits = 2), sep = ""),
                              paste( "p-value = ",      format(test2$p.value,   digits = 2), sep = ""),
			      paste( "p-value = ",      format(test3$p.value,   digits = 2), sep = "")  ),
				col=c("blue", "black", "green"), pch=16 )

    dev.off()
}

discrete <- function (ms, scrm, scrm_p10000, scrm_p50000, title_i){
    test1 = ks.test( ms, scrm )
    test2 = ks.test( ms, scrm_p50000)
    test3 = ks.test( ms, scrm_p10000)
    xrange = range ( c(ms, scrm, scrm_p10000, scrm_p50000) )
    pdf ( paste( title_i,".pdf", sep = "" ) )
    plot( ecdf(ms), xlim = xrange, col="red", main = paste( title_i ) , xlab = title_i)
    plot( ecdf(scrm), add=TRUE, lty="dashed", col="blue")
    plot( ecdf(scrm_p10000), add=TRUE, lty="dashed", col="green")
    plot( ecdf(scrm_p50000), add=TRUE, lty="dashed", col="black")
    legend( "topleft" , c( "window 10000", "window 50000", "scrm", "ms"), col=c( "green", "black",  "blue", "red" ), pch=16)
    legend( "bottomright", c( paste( "p-value = ",      format(test1$p.value,   digits = 2), sep = ""),
                              paste( "p-value = ",      format(test2$p.value,   digits = 2), sep = ""),
			      paste( "p-value = ",      format(test3$p.value,   digits = 2), sep = "")  ),
				col=c("blue", "black", "green"), pch=16 )

    dev.off()
}

#ms     = read.table( "ms_stat"   ) 
#scrm   = read.table( "scrm_stat" ) 
#scrm_p = read.table( "scrmprune_stat") 

ms_     = read.table( "ms_stat")
scrm_   = read.table( "scrm_stat")
scrm_p0_ = read.table( "scrmprune0_stat")
#scrm_=scrm_p0_
scrm_p10000_ = read.table( "scrmprune10000_stat")
scrm_p50000_ = read.table( "scrmprune50000_stat")

discrete ( ms_$V2, scrm_$V2, scrm_p10000_$V2, scrm_p50000_$V2, "Pairwise_difference")
discrete ( ms_$V4, scrm_$V4, scrm_p10000_$V4, scrm_p50000_$V4, "ss")

#continuous ( ms_$V2, scrm_$V2, scrm_p10000_$V2, scrm_p50000_$V2, "Pairwise_difference")
#continuous ( ms_$V4, scrm_$V4, scrm_p10000_$V4, scrm_p50000_$V4, "ss")
continuous ( ms_$V6, scrm_$V6, scrm_p10000_$V6, scrm_p50000_$V6, "TajimaD")
continuous ( ms_$V8, scrm_$V8, scrm_p10000_$V8, scrm_p50000_$V8, "ThetaH")
continuous ( ms_$V10, scrm_$V10, scrm_p10000_$V10, scrm_p50000_$V10, "H")

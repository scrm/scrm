rm(list=ls())

tmrca_moment <- function ( ms, scrm, scrm_p){
    ms_scrm_test = list ( ks.test( ms$V1, scrm$V1 ),
                          ks.test( ms$V2, scrm$V2 ),
                          ks.test( ms$V3, scrm$V3 ),
                          ks.test( ms$V4, scrm$V4 ) )
    
    ms_scrmp_test = list ( ks.test( ms$V1, scrm_p$V1 ),
                           ks.test( ms$V2, scrm_p$V2 ),
                           ks.test( ms$V3, scrm_p$V3 ),
                           ks.test( ms$V4, scrm_p$V4 ) )    
    
    titles = c("1st", "2nd", "3rd", "4th")
    
    pdf ( "TMRCA_moment_KStest.pdf" )
    par(mfrow=c(2,2))
    for ( i in c(1:4) ) {
        ms_i    = ms[,i]
        scrm_i  = scrm[,i]
        scrm_pi = scrm_p[,i]
        xrange = range ( c(ms_i, scrm_i, scrm_pi) )
        plot( ecdf(ms_i), xlim = xrange, col="red", main = paste( titles[i], "Moment" ) )
        plot( ecdf(scrm_i), add=TRUE, lty="dashed", col="blue")
        plot( ecdf(scrm_pi), add=TRUE, lty="dashed", col="green")
        legend( "bottomright", c( paste( "Test-stats 1 = ", format(ms_scrm_test[[i]]$statistic, digits = 2), sep = ""), 
                                  paste( "p-value = ",      format(ms_scrm_test[[i]]$p.value,   digits = 2), sep = ""),
                                  paste( "Test-stats 2 = ", format(ms_scrmp_test[[i]]$statistic,digits = 2), sep = ""), 
                                  paste( "p-value = ",      format(ms_scrmp_test[[i]]$p.value,  digits = 2), sep = "")  ) )
        legend( "topleft" , c( "scrm 50000 exact", "scrm", "ms"), col=c( "green", "blue", "red" ), pch=16)
    }
    dev.off()
}

ms_     = read.table( "ms_moment" ,         sep = ",") 
scrm_   = read.table( "scrm_moment" ,       sep = ",") 
scrm_p_ = read.table( "scrmprune_moment" , sep = ",") 

tmrca_moment ( ms_, scrm_, scrm_p_)

#cat(paste(${nsam},${r},\"|\",
#format(m1test\$statistic,digits=4),format(m1test\$p.value,scientific = TRUE),\"|\",
#format(m2test\$statistic,digits=4),format(m2test\$p.value,scientific = TRUE),\"|\",
#format(m3test\$statistic,digits=4),format(m3test\$p.value,scientific = TRUE),\"|\",
#format(m4test\$statistic,digits=4),format(m4test\$p.value,scientific = TRUE),
#sep=\"\t\"),file=\"${compareTMRCA}\",append=TRUE);cat(\"\n\",file=\"${compareTMRCA}\",append=TRUE);

rm( list=ls() )

continuous <- function (ms, scrm, scrm_p, title_i){
    test1 = ks.test( ms, scrm )
    test2 = ks.test( ms, scrm_p)
    xrange = range ( c(ms, scrm, scrm_p) )
    pdf ( paste( title_i,".pdf", sep = "" ) )
    plot( ecdf(ms), xlim = xrange, col="red", main = paste( title_i ) , xlab = title_i)
    plot( ecdf(scrm), add=TRUE, lty="dashed", col="blue")
    plot( ecdf(scrm_p), add=TRUE, lty="dashed", col="green")
    legend( "topleft" , c( "scrm 50000 exact", "scrm", "ms"), col=c( "green", "blue", "red" ), pch=16)
    legend( "bottomright", c( paste( "Test-stats 1 = ", format(test1$statistic, digits = 2), sep = ""), 
                              paste( "p-value = ",      format(test1$p.value,   digits = 2), sep = ""),
                              paste( "Test-stats 2 = ", format(test2$statistic,digits = 2), sep = ""), 
                              paste( "p-value = ",      format(test2$p.value,  digits = 2), sep = "")  ) )

    dev.off()
}
#echo "rm(list=ls());
##source(\"../fun_src.r\");
#figuretitle=scan(\"figuretitle\",what=\"\");
#currentcase=scan(\"current_case\",what=\"\");
#msdata=read.table(paste(\"ms\",\"data\",sep=\"\"))\$V1;
#scrmdata=read.table(paste(\"scrm\",\"data\",sep=\"\"))\$V1;
##ee=1#ee_tmrca(${nsam});
##sdv=1#sd_tmrca(${nsam});
#test=ks.test(msdata,scrmdata)
#pdf(paste(currentcase,figuretitle,\".pdf\",sep=\"\"));
#plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=currentcase)
#plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
#legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,digits=4),sep=\"\")))
#legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
#dev.off();
##cat(paste(currentcase,\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
#cat(paste(currentcase, figuretitle , \"\n\",\"|\",
#format(mean(msdata),scientific = TRUE),format(sd(msdata),scientific = TRUE),\"||\",
#format(mean(scrmdata),scientific = TRUE),format(sd(scrmdata),scientific = TRUE),\"|\",format(test\$statistic,scientific = TRUE),format(test\$p.value,scientific = TRUE), 
#sep=\"\t\"),file=\"${compareMIG}\",append=TRUE);cat(\"\n\",file=\"${compareMIG}\",append=TRUE);" > ks.r


discrete <- function (ms, scrm, scrm_p, title_i){
    xrange = range ( c(ms, scrm, scrm_p) )
#    mstable    = table( factor (ms ,    levels=xrange[1]:xrange[2]))
##    print(mstable)
#    scrmtable  = table( factor (scrm,   levels=xrange[1]:xrange[2]))
#    scrmptable = table( factor (scrm_p, levels=xrange[1]:xrange[2]))
    mstable    = table(ms)
    print(mstable)
    scrmtable  = table( scrm)
    scrmptable = table( scrm_p )

    yrange = range ( c(mstable, scrmtable, scrmptable) )
    pdf ( paste( title_i,".pdf", sep = "" ) )
    plot ( xrange, yrange, type = "n" , ylab = "Frequency", xlab = title_i)
    points ( as.numeric( names(mstable  )  ),  mstable   ,   pch=1, col = "red" )
    points ( as.numeric( names(scrmtable)  ),  scrmtable ,   pch=1, col = "blue"  )
    points ( as.numeric( names(scrmptable) ),  scrmptable ,  pch=1, col = "green" )
    legend( "topleft" , c( "scrm 50000 exact", "scrm", "ms"), col=c( "green", "blue", "red" ), pch=16)
    dev.off()
}

ms     = read.table( "ms_stat"   ) 
scrm   = read.table( "scrm_stat" ) 
scrm_p = read.table( "scrmprune_stat") 


#discrete ( ms$V2, scrm$V2, scrm_p$V2, "Pairwise_difference")
#discrete ( ms$V4, scrm$V4, scrm_p$V4, "ss")

continuous ( ms$V2, scrm$V2, scrm_p$V2, "Pairwise_difference")
continuous ( ms$V4, scrm$V4, scrm_p$V4, "ss")
continuous ( ms$V6, scrm$V6, scrm_p$V6, "TajimaD")
continuous ( ms$V8, scrm$V8, scrm_p$V8, "ThetaH")
continuous ( ms$V10, scrm$V10, scrm_p$V10, "H")
#figuretitle=scan(\"figuretitle\",what=\"\");
#currentcase=scan(\"current_case\",what=\"\");
#msdata=read.table(paste(\"ms\",\"data\",sep=\"\"))\$V1;
#scrmdata=read.table(paste(\"scrm\",\"data\",sep=\"\"))\$V1;
#		#ee=ee_seg(${nsam},${t});
#		#sdv=sd_seg_norecomb(${nsam},${t});
#		datamax=max(msdata,scrmdata);
#		mstable=table(factor(msdata,levels=1:datamax))
#		scrmtable=table(factor(scrmdata,levels=1:datamax))
#pdf(paste(currentcase,figuretitle,\".pdf\",sep=\"\"));
#plot(as.numeric(names(mstable)), mstable/length(msdata),pch=16,col=\"red\",ylab=\"Frequency\",xlab=figuretitle);
#points(as.numeric(names(scrmtable)), scrmtable/length(scrmdata),pch=16,col=\"blue\")

#ms_newtable=table(msdata);
#scrm_mstable=table(factor(scrmdata,levels=names(table(msdata))));
#combined_scrm_ms_test=chisq.test(cbind(scrm_mstable, ms_newtable));

#scrm_newtable=table(scrmdata);
#ms_scrmtable=table(factor(msdata,levels=names(table(scrmdata))));
#combined_ms_scrm_test=chisq.test(cbind(scrm_newtable, ms_scrmtable));

#legend(\"topright\",c(paste(\"Test 1 Statistics = \",combined_scrm_ms_test\$statistic,sep=\"\"), paste(\"p-value = \",format(combined_scrm_ms_test\$p.value,digits=4),sep=\"\"),paste(\"Test 2 Statistics = \",combined_ms_scrm_test\$statistic,sep=\"\"), paste(\"p-value = \",format(combined_ms_scrm_test\$p.value,digits=4),sep=\"\")))
#legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)

#dev.off();

##cat(paste(${nsam},${t},\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
#cat(paste(currentcase, figuretitle , \"\n\",\"|\",
#format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
#format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"||\",
#format(combined_scrm_ms_test\$statistic,digits=4),format(combined_scrm_ms_test\$p.value,scientific = TRUE),\"||\",
#format(combined_ms_scrm_test\$statistic,digits=4),format(combined_ms_scrm_test\$p.value,scientific = TRUE),
#sep=\"\t\"),file=\"${compareMIG}\",append=TRUE);cat(\"\n\",file=\"${compareMIG}\",append=TRUE);
#" > chisq.r

#echo "rm(list=ls());
##source(\"../fun_src.r\");
#figuretitle=scan(\"figuretitle\",what=\"\");
#currentcase=scan(\"current_case\",what=\"\");
#msdata=read.table(paste(\"ms\",\"data\",sep=\"\"))\$V1;
#scrmdata=read.table(paste(\"scrm\",\"data\",sep=\"\"))\$V1;
##ee=1#ee_tmrca(${nsam});
##sdv=1#sd_tmrca(${nsam});
#test=ks.test(msdata,scrmdata)
#pdf(paste(currentcase,figuretitle,\".pdf\",sep=\"\"));
#plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=currentcase)
#plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
#legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,digits=4),sep=\"\")))
#legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
#dev.off();
##cat(paste(currentcase,\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
#cat(paste(currentcase, figuretitle , \"\n\",\"|\",
#format(mean(msdata),scientific = TRUE),format(sd(msdata),scientific = TRUE),\"||\",
#format(mean(scrmdata),scientific = TRUE),format(sd(scrmdata),scientific = TRUE),\"|\",format(test\$statistic,scientific = TRUE),format(test\$p.value,scientific = TRUE), 
#sep=\"\t\"),file=\"${compareMIG}\",append=TRUE);cat(\"\n\",file=\"${compareMIG}\",append=TRUE);" > ks.r

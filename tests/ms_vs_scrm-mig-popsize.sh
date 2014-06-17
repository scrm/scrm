#!/bin/bash

#Compare summary statistics of ms and scrm for migration

mkdir test-mig-popsize
cd test-mig-popsize
rm -rf *pdf


rep=10000

## compare TMRCA
compareMIG=compareMIG
rm ${comparePop}


#rm ${compareTMRCA}
#rm ${compareBL}
#echo -e "compare TMRCA for ${rep} replicates 
#\t|\t ms \t\t|\t scrm\t\t|\tKS test
#Case\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value" >${compareTMRCA}

#echo -e "compare BL for ${rep} replicates 
#\t|\t ms \t\t|\t scrm\t\t|\tKS test
#Case\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value" >${compareBL}

theta=10

		echo "rm(list=ls());
		#source(\"../fun_src.r\");
figuretitle=scan(\"figuretitle\",what=\"\");
currentcase=scan(\"current_case\",what=\"\");
msdata=read.table(paste(\"ms\",\"data\",sep=\"\"))\$V1;
scrmdata=read.table(paste(\"scrm\",\"data\",sep=\"\"))\$V1;
		#ee=ee_seg(${nsam},${t});
		#sdv=sd_seg_norecomb(${nsam},${t});
		datamax=max(msdata,scrmdata);
		mstable=table(factor(msdata,levels=1:datamax))
		scrmtable=table(factor(scrmdata,levels=1:datamax))
pdf(paste(currentcase,figuretitle,\".pdf\",sep=\"\"));
plot(as.numeric(names(mstable)), mstable/length(msdata),pch=16,col=\"red\",ylab=\"Frequency\",xlab=figuretitle);
points(as.numeric(names(scrmtable)), scrmtable/length(scrmdata),pch=16,col=\"blue\")

ms_newtable=table(msdata);
scrm_mstable=table(factor(scrmdata,levels=names(table(msdata))));
combined_scrm_ms_test=chisq.test(cbind(scrm_mstable, ms_newtable));

scrm_newtable=table(scrmdata);
ms_scrmtable=table(factor(msdata,levels=names(table(scrmdata))));
combined_ms_scrm_test=chisq.test(cbind(scrm_newtable, ms_scrmtable));

legend(\"topright\",c(paste(\"Test 1 Statistics = \",combined_scrm_ms_test\$statistic,sep=\"\"), paste(\"p-value = \",format(combined_scrm_ms_test\$p.value,digits=4),sep=\"\"),paste(\"Test 2 Statistics = \",combined_ms_scrm_test\$statistic,sep=\"\"), paste(\"p-value = \",format(combined_ms_scrm_test\$p.value,digits=4),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)

dev.off();

#cat(paste(${nsam},${t},\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
cat(paste(currentcase, figuretitle , \"\n\",\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"||\",
format(combined_scrm_ms_test\$statistic,digits=4),format(combined_scrm_ms_test\$p.value,scientific = TRUE),\"||\",
format(combined_ms_scrm_test\$statistic,digits=4),format(combined_ms_scrm_test\$p.value,scientific = TRUE),
sep=\"\t\"),file=\"${compareMIG}\",append=TRUE);cat(\"\n\",file=\"${compareMIG}\",append=TRUE);
" > chisq.r

echo "rm(list=ls());
#source(\"../fun_src.r\");
figuretitle=scan(\"figuretitle\",what=\"\");
currentcase=scan(\"current_case\",what=\"\");
msdata=read.table(paste(\"ms\",\"data\",sep=\"\"))\$V1;
scrmdata=read.table(paste(\"scrm\",\"data\",sep=\"\"))\$V1;
#ee=1#ee_tmrca(${nsam});
#sdv=1#sd_tmrca(${nsam});
test=ks.test(msdata,scrmdata)
pdf(paste(currentcase,figuretitle,\".pdf\",sep=\"\"));
plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=currentcase)
plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,digits=4),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
dev.off();
#cat(paste(currentcase,\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
cat(paste(currentcase, figuretitle , \"\n\",\"|\",
format(mean(msdata),scientific = TRUE),format(sd(msdata),scientific = TRUE),\"||\",
format(mean(scrmdata),scientific = TRUE),format(sd(scrmdata),scientific = TRUE),\"|\",format(test\$statistic,scientific = TRUE),format(test\$p.value,scientific = TRUE), 
sep=\"\t\"),file=\"${compareMIG}\",append=TRUE);cat(\"\n\",file=\"${compareMIG}\",append=TRUE);" > ks.r

echo "rm(list=ls());
#source(\"../fun_src.r\");
figuretitle=scan(\"figuretitle\",what=\"\");
currentcase=scan(\"current_case\",what=\"\");
msdata=read.table(\"mstmrca\")\$V1;
scrmdata=read.table(\"scrmtmrca\")\$V1;
#ee=ee_tmrca(${nsam});
#sdv=sd_tmrca(${nsam});
test=ks.test(msdata,scrmdata)
pdf(paste(currentcase,figuretitle,\".pdf\",sep=\"\"));
plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=currentcase)
plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,digits=4),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
dev.off();
cat(paste(currentcase,figuretitle , \"\n\",\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"|\",test\$statistic,format(test\$p.value,digits=4), 
sep=\"\t\"),file=\"${compareMIG}\",append=TRUE);cat(\"\n\",file=\"${compareMIG}\",append=TRUE);" > tmrca.r

echo "rm(list=ls());
#source(\"../fun_src.r\");
figuretitle=scan(\"figuretitle\",what=\"\");
currentcase=scan(\"current_case\",what=\"\");
msdata=read.table(\"msbl\")\$V1;
scrmdata=read.table(\"scrmbl\")\$V1;
test=ks.test(msdata,scrmdata)
pdf(paste(currentcase,figuretitle,\".pdf\",sep=\"\"));
plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=currentcase)
plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,digits=4),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
dev.off();
cat(paste(currentcase,figuretitle , \"\n\",\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"|\",test\$statistic,format(test\$p.value,digits=4), 
sep=\"\t\"),file=\"${compareMIG}\",append=TRUE);cat(\"\n\",file=\"${compareMIG}\",append=TRUE);" > bl.r


foo(){
    cut -f 2 mstime > mstmrca
	cut -f 2 scrmtime > scrmtmrca
	echo "TMRCA" > figuretitle
	R CMD BATCH tmrca.r

    cut -f 3 mstime > msbl
	cut -f 3 scrmtime > scrmbl
	echo "BL" > figuretitle
	R CMD BATCH bl.r

	cut -f 6 ms_stats > msdata
	cut -f 6 scrm_stats > scrmdata
	echo "Tajima_D" > figuretitle
	R CMD BATCH ks.r

	cut -f 2 ms_stats > msdata
	cut -f 2 scrm_stats > scrmdata
	echo "Pairewise_difference" > figuretitle
	R CMD BATCH chisq.r

	cut -f 8 ms_stats > msdata
	cut -f 8 scrm_stats > scrmdata
	echo "theta_H" > figuretitle
	R CMD BATCH chisq.r

	cut -f 10 ms_stats > msdata
	cut -f 10 scrm_stats > scrmdata
	echo "H" > figuretitle
	R CMD BATCH chisq.r
	}
	
	
	
#case 1 
#2 sub population, 2 samples from each subpopulation, mutation rate is 5
echo "2groups2sam2sam_mig5" > current_case
rm ms* scrm*

#ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}

ms 4 ${rep} -t ${theta} -I 2 2 2 5.0 -eN 0.4 10.01 -eN 1 0.01 -T -L > msout
scrm 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -eN 0.4 10.01 -eN 1 0.01 -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo

#case 2 
#2 sub population, 6 samples from subpopulation 1, and 5 samples from subpopulation 2, with rate 10 from 1 to 2, and rate 5 from 2 to 1

echo "2groups6sam5sam_mig_x_10_5_x" > current_case
rm ms* scrm*
	
ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 10.0 5.0 x -eN .5 0.01 -T -L > msout
scrm 4 ${rep} -t ${theta} -I 2 2 2 -ma x 10.0 5.0 x -eN .5 0.01 -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo


#case 3
#3 sub population, 10 samples from subpopulation 1, and 4 samples from subpopulation 2, and 1 sample from the third with rate 5

echo "3groups10sam4sam1sam_mig_sym10" > current_case
rm ms* scrm*
	
ms 15 ${rep} -t ${theta} -I 3 10 4 1 10.0 -eN 0.5 10.0 -T -L > msout
scrm 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 5.0 5.0 5.0 x 5.0 5.0 5.0 x -eN 0.5 10.0 -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo

#case 4
#3 sub population, 10 samples from subpopulation 1, and 4 samples from subpopulation 2, and 1 sample from the third with rate 5

echo "3groups10sam4sam1sam_mig_offdiag5" > current_case
rm ms* scrm*

	
ms 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 5.0 5.0 5.0 x 5.0 5.0 5.0 x -eN 0.8 15  -T -L > msout
scrm 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 5.0 5.0 5.0 x 5.0 5.0 5.0 x -eN 0.8 15  -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo



#case 5
#3 sub population, 10 samples from subpopulation 1, and 4 samples from subpopulation 2, and 1 sample from the third with rate matrix on manual page 5

echo "3groups10sam4sam1sam_mig_x123x456x" > current_case

rm ms* scrm*
	
ms 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x -eN 1 .1 -eN 3 10 -T -L > msout
scrm 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x -eN 1 .1 -eN 3 10 -T -L > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep "time:" > mstime

cat scrmout | sample_stats > scrm_stats
cat scrmout | grep "time:" >  scrmtime

foo

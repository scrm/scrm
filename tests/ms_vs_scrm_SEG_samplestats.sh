#!/bin/bash

mkdir test-SEG
cd test-SEG
rm *pdf


rep=100000

## compare SEG data
compareSEG=compareSEG
rm ${compareSEG}

theta=10


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
format(mean(msdata),scientific = TRUE),format(sd(msdata),scientific = TRUE),\"|\",
format(mean(scrmdata),scientific = TRUE),format(sd(scrmdata),scientific = TRUE),\"|\",format(test\$statistic,scientific = TRUE),format(test\$p.value,scientific = TRUE), 
sep=\"\t\"),file=\"${compareSEG}\",append=TRUE);cat(\"\n\",file=\"${compareSEG}\",append=TRUE);" > dummy.r

foo(){
	cut -f 2 ms_stats > msdata
	cut -f 2 scrm_stats > scrmdata
	echo "Pairewise_difference" > figuretitle
	R CMD BATCH dummy.r
	
	cut -f 6 ms_stats > msdata
	cut -f 6 scrm_stats > scrmdata
	echo "Tajima_D" > figuretitle
	R CMD BATCH dummy.r
	
	cut -f 8 ms_stats > msdata
	cut -f 8 scrm_stats > scrmdata
	echo "theta_H" > figuretitle
	R CMD BATCH dummy.r
	
	cut -f 10 ms_stats > msdata
	cut -f 10 scrm_stats > scrmdata
	echo "H" > figuretitle
	R CMD BATCH dummy.r
	}

#case 1 
echo "4_samples" > current_case

ms 4 ${rep} -t ${theta} | sample_stats > ms_stats
scrm 4 ${rep} -t ${theta} | sample_stats > scrm_stats

foo


#case 2

echo "5_samples" > current_case

ms 5 ${rep} -t ${theta} | sample_stats > ms_stats
scrm 5 ${rep} -t ${theta} | sample_stats > scrm_stats

foo

#case 3

#2 sub population, 6 samples from subpopulation 1, and 5 samples from subpopulation 2, with rate 10 from 1 to 2, and rate 5 from 2 to 1

prefix="2groups6sam5sam_mig_x_10_5_x"
echo ${prefix} > current_case
	
ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 10.0 5.0 x | sample_stats > ms_stats
scrm 4 ${rep} -t ${theta} -I 2 2 2 -ma x 10.0 5.0 x | sample_stats > scrm_stats

foo

#case 4
#3 sub population, 10 samples from subpopulation 1, and 4 samples from subpopulation 2, and 1 sample from the third with rate matrix on manual page 5

prefix="3groups10sam4sam1sam_mig_x123x456x"
echo ${prefix} > current_case
	
ms 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x | sample_stats > ms_stats
scrm 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x | sample_stats > scrm_stats

foo

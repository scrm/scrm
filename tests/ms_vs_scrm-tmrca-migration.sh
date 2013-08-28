#!/bin/bash

#Compare summary statistics of ms and scrm for TMRCA for migration

mkdir test-tmrca-mig
cd test-tmrca-mig
rm *pdf


rep=100000

## compare TMRCA
compareTMRCA=compareTMRCA
compareBL=compareBL

rm ${compareTMRCA}
rm ${compareBL}
echo -e "compare TMRCA for ${rep} replicates 
\t|\t ms \t\t|\t scrm\t\t|\tKS test
Case\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value" >${compareTMRCA}

echo -e "compare BL for ${rep} replicates 
\t|\t ms \t\t|\t scrm\t\t|\tKS test
Case\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value" >${compareBL}

theta=10

echo "rm(list=ls());
#source(\"../fun_src.r\");
currentcase=scan(\"current_case\",what=\"\");
msdata=read.table(paste(\"ms\",currentcase,\"tmrca\",sep=\"\"))\$V1;
scrmdata=read.table(paste(\"scrm\",currentcase,\"tmrca\",sep=\"\"))\$V1;
#ee=1#ee_tmrca(${nsam});
#sdv=1#sd_tmrca(${nsam});
test=ks.test(msdata,scrmdata)
pdf(paste(currentcase,\"TMRCA-KStest.pdf\",sep=\"\"));
plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=currentcase)
plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,digits=4),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
dev.off();
#cat(paste(currentcase,\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
cat(paste(currentcase,\"\n\",\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"|\",test\$statistic,format(test\$p.value,digits=4), 
sep=\"\t\"),file=\"${compareTMRCA}\",append=TRUE);cat(\"\n\",file=\"${compareTMRCA}\",append=TRUE);
rm(list=ls());
#source(\"../fun_src.r\");
currentcase=scan(\"current_case\",what=\"\");
msdata=read.table(paste(\"ms\",currentcase,\"bl\",sep=\"\"))\$V1;
scrmdata=read.table(paste(\"scrm\",currentcase,\"bl\",sep=\"\"))\$V1;
#ee=1#ee_tmrca(${nsam});
#sdv=1#sd_tmrca(${nsam});
test=ks.test(msdata,scrmdata)
pdf(paste(currentcase,\"BL-KStest.pdf\",sep=\"\"));
plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=currentcase)
plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,digits=4),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
dev.off();
#cat(paste(currentcase,\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
cat(paste(currentcase,\"\n\",\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"|\",test\$statistic,format(test\$p.value,digits=4), 
sep=\"\t\"),file=\"${compareBL}\",append=TRUE);cat(\"\n\",file=\"${compareBL}\",append=TRUE);" > dummy.r

#case 1 
#2 sub population, 2 samples from each subpopulation, mutation rate is 5
prefix="2groups2sam2sam_mig5"
Trees=${prefix}Trees
bl=${prefix}bl
tmrca=${prefix}tmrca
	
echo ${prefix} > current_case
	
#ms 4 ${rep} -t ${theta} -I 2 2 2 5.0 -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}
ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}
hybrid-Lambda -gt ms${Trees} -tmrca ms${tmrca}
hybrid-Lambda -gt ms${Trees} -bl ms${bl}

scrm 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > scrm${Trees}
hybrid-Lambda -gt scrm${Trees} -tmrca scrm${tmrca}
hybrid-Lambda -gt scrm${Trees} -bl scrm${bl}

R CMD BATCH dummy.r
rm ms* scrm*

#case 2 
#2 sub population, 6 samples from subpopulation 1, and 5 samples from subpopulation 2, with rate 10 from 1 to 2, and rate 5 from 2 to 1

prefix="2groups6sam5sam_mig_x_10_5_x"
Trees=${prefix}Trees
bl=${prefix}bl
tmrca=${prefix}tmrca
	
echo ${prefix} > current_case
	
ms 4 ${rep} -t ${theta} -I 2 2 2 5.0 -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}
hybrid-Lambda -gt ms${Trees} -tmrca ms${tmrca}
hybrid-Lambda -gt ms${Trees} -bl ms${bl}

scrm 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > scrm${Trees}
hybrid-Lambda -gt scrm${Trees} -tmrca scrm${tmrca}
hybrid-Lambda -gt scrm${Trees} -bl scrm${bl}

R CMD BATCH dummy.r
rm ms* scrm*

#case 3 and 4 are from ms
#case 3
#3 sub population, 10 samples from subpopulation 1, and 4 samples from subpopulation 2, and 1 sample from the third with rate 5

prefix="3groups10sam4sam1sam_mig_5"
Trees=${prefix}Trees
bl=${prefix}bl
tmrca=${prefix}tmrca
	
echo ${prefix} > current_case
	
ms 15 ${rep} -t ${theta} -I 3 10 4 1 5.0 -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}
hybrid-Lambda -gt ms${Trees} -tmrca ms${tmrca}
hybrid-Lambda -gt ms${Trees} -bl ms${bl}

scrm 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 5.0 5.0 5.0 x 5.0 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > scrm${Trees}
hybrid-Lambda -gt scrm${Trees} -tmrca scrm${tmrca}
hybrid-Lambda -gt scrm${Trees} -bl scrm${bl}

R CMD BATCH dummy.r
rm ms* scrm*


#case 4
#3 sub population, 10 samples from subpopulation 1, and 4 samples from subpopulation 2, and 1 sample from the third with rate matrix on manual page 5

prefix="3groups10sam4sam1sam_mig_x123x456x"
Trees=${prefix}Trees
bl=${prefix}bl
tmrca=${prefix}tmrca
	
echo ${prefix} > current_case
	
ms 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}
hybrid-Lambda -gt ms${Trees} -tmrca ms${tmrca}
hybrid-Lambda -gt ms${Trees} -bl ms${bl}

scrm 15 ${rep} -t ${theta} -I 3 10 4 1 -ma x 1.0 2.0 3.0 x 4.0 5.0 6.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > scrm${Trees}
hybrid-Lambda -gt scrm${Trees} -tmrca scrm${tmrca}
hybrid-Lambda -gt scrm${Trees} -bl scrm${bl}

R CMD BATCH dummy.r
rm ms* scrm*

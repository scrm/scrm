#!/bin/bash

#Compare summary statistics of ms and scrm for TMRCA, number of mutation, and number of recombination

msNsample=(3 4 5)

rep=10000

## compare TMRCA
compareTMRCA=compareTMRCA
#rm ${compareTMRCA}
#theta=10
npop=1000000
for nsam in "${msNsample[@]}"
	do
	prefix=${nsam}sample
	out=${prefix}out
	freq=${prefix}freq
	Trees=${prefix}Trees
	ms ${nsam} ${rep} -T | tail -n +4 | grep -v "//" > ms${out}
	cat ms${out} | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}
	hybrid-Lambda -gt ms${Trees} -fF ms${freq}
	scrm ${nsam} ${rep} -T | tail -n +4 | grep -v "//" > scrm${out}
	cat scrm${out} | grep ";" | sed -e 's/\[.*\]//g' > scrm${Trees}
	hybrid-Lambda -gt scrm${Trees} -fF scrm${freq}
	
	#echo "rm(list=ls());msdata=read.table(\"ms${freq}\")\$V1;scrmdata=read.table(\"scrm${freq}\")\$V1;cat(paste(${nsam},mean(msdata),sd(msdata),sd(msdata)/sqrt(length(msdata)),mean(scrmdata),sd(scrmdata),sd(scrmdata)/sqrt(length(scrmdata)),sep=\"\t\"),file=\"${compareTMRCA}\",append=TRUE);cat(\"\n\",file=\"${compareTMRCA}\",append=TRUE);" > dummy.r
	#R CMD BATCH dummy.r
	#rm ms${out} ms${Trees} ms${freq} scrm${out} scrm${Trees} scrm${freq} 
	done


#!/bin/bash

###########################

mst=(10 20 50 100)
msNsample=(2 3 7 10)
rep=100
#npop=20000

## compare number of segregating sites
compareSEG=compareSEG
rm ${compareSEG}
for t in "${mst[@]}"
	do
	for nsam in "${msNsample[@]}"
		do
		prefix=${nsam}sample${t}mut
		out=${prefix}out
		nseg=${prefix}NumOfSeg
		ms ${nsam} ${rep} -t ${t} -T | tail -n +4 | grep -v "//" > ms${out}
		cat ms${out} | grep "segsites" | sed -e "s/segsites: //" > ms${nseg}
	
		scrm ${nsam} ${rep} -t ${t} | tail -n +4 | grep -v "//" > scrm${out}
		cat scrm${out} | grep "segsites" | sed -e "s/segsites: //" > scrm${nseg}
		
		echo "rm(list=ls());msdata=read.table(\"ms${nseg}\")\$V1;scrmdata=read.table(\"scrm${nseg}\")\$V1;cat(paste(${nsam},${t},mean(msdata),sd(msdata),sd(msdata)/sqrt(length(msdata)),mean(scrmdata),sd(scrmdata),sd(scrmdata)/sqrt(length(scrmdata)),sep=\"\t\"),file=\"${compareSEG}\",append=TRUE);cat(\"\n\",file=\"${compareSEG}\",append=TRUE);" > dummy.r
		R CMD BATCH dummy.r
		rm ms${out} ms${nseg} scrm${out} scrm${nseg}
		done
	done


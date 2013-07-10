#!/bin/bash

###########################

mst=(10 20 50 100)
msNsample=(2 3 7 10)
rep=1000
#npop=20000

## compare number of segregating sites
compareSEG=compareSEG
echo -e "compare number of mutations for ${rep} replicates \n\t\t|\t\t\t|\t\t ms \t\t|\t\t scrm\nNsam\ttheta\t|\tmean\tstdv\t|\tmean\tstdv\tstd err\t|\tmean\tstdv \tstd err" >${compareSEG}

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
		
		echo "rm(list=ls());
		source(\"fun_src.r\");
		msdata=read.table(\"ms${nseg}\")\$V1;
		scrmdata=read.table(\"scrm${nseg}\")\$V1;
		ee=ee_seg(${nsam},${t});
		sdv=sd_seg_norecomb(${nsam},${t});
		cat(paste(${nsam},${t},\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),format(sd(msdata)/sqrt(length(msdata)),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),format(sd(scrmdata)/sqrt(length(scrmdata)),digits=4),
sep=\"\t\"),file=\"${compareSEG}\",append=TRUE);cat(\"\n\",file=\"${compareSEG}\",append=TRUE);" > dummy.r
		R CMD BATCH dummy.r
		#rm ms${out} ms${nseg} scrm${out} scrm${nseg}
		done
	done


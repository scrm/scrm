#!/bin/bash

## compare times of genealogy changes
msr=(10 20 10 50 100)
mst=(10 20 50 100 10)

rep=100
seqlen=100000
msNsample=(2 3 7 10)
compareRECOMB=compareSEG-RECOMB
echo -e "compare number of mutations for ${rep} replicates \n\t\t\t|\t\t\t|\t\t ms \t\t|\t\t scrm\nNsam\ttheta\trho\t|\tmean\tstdv\t|\tmean\tstdv\tstd err\t|\tmean\tstdv \tstd err" > ${compareRECOMB}

for t in "${mst[@]}"
	do
	for r in "${msr[@]}"
		do
		for nsam in "${msNsample[@]}"
			do
			prefix=${nsam}sample${r}rho${t}theta
			out=${prefix}out
			segrecomb=${prefix}segRecomb

			ms ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} -T | tail -n +4 | grep -v "//" > ms${out}
			cat ms${out} | grep "segsites" | sed -e "s/segsites: //" > ms${segrecomb}
	
			scrm ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} | tail -n +4 | grep -v "//" > scrm${out}
			cat scrm${out} | grep "segsites" | sed -e "s/segsites: //" > scrm${segrecomb}

			
			echo "rm(list=ls());
			source(\"fun_src.r\");
			msdata=read.table(\"ms${segrecomb}\")\$V1;
			scrmdata=read.table(\"scrm${segrecomb}\")\$V1;
			ee=ee_seg(${nsam},${r});
			sdv=sd_seg_recomb(${nsam},${t},${r});
			cat(paste(${nsam},${t},${r},\"|\",
format(ee,digits=4),format(sdv,digits=4),\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),format(sd(msdata)/sqrt(length(msdata)),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),format(sd(scrmdata)/sqrt(length(scrmdata)),digits=4),
sep=\"\t\"),file=\"${compareRECOMB}\",append=TRUE);cat(\"\n\",file=\"${compareRECOMB}\",append=TRUE);" > dummy.r
			R CMD BATCH dummy.r
			rm ms${out} ms${segrecomb} scrm${out} scrm${segrecomb} 
			done
		done
	done

	
	

#!/bin/bash

## compare times of genealogy changes
msr=(10 20 10 50)
rep=1000
seqlen=100000
msNsample=(2 3 7 10 20)
compareRECOMB=compareRECOMB
echo -e "compare number of recombination for ${rep} replicates \n\t\t|\t\t\t|\t\t ms \t\t|\t\t scrm\nNsam\trho\t|\tmean\tstdv\t|\tmean\tstdv\tstd err\t|\tmean\tstdv \tstd err" >${compareRECOMB}
#for t in "${mst[@]}"
#	do
	for r in "${msr[@]}"
		do
		for nsam in "${msNsample[@]}"
			do
			prefix=${nsam}sample${r}recomb
			out=${prefix}out
			recomb=${prefix}Recomb
			ms ${nsam} ${rep} -r ${r} ${seqlen} -T | tail -n +4  > ms${out}
			rm xx*
			csplit -n 1 ms${out} '/^///' "{*}"
			rm ms${recomb}
			for file in $(seq 1 1 ${rep})
				do 
				cat xx${file} | grep ";" | wc -l >> ms${recomb}
				done
				
				
			scrm ${nsam} ${rep} -r ${r} ${seqlen} -T | tail -n +4  > scrm${out}
			
			rm xx*
			csplit -n 1 scrm${out}  '/^///' "{*}"
			rm scrm${recomb}
			for file in $(seq 1 1 ${rep})
				do 
				cat xx${file} | grep ";" | wc -l >> scrm${recomb}
				done
			echo "rm(list=ls());
			source(\"fun_src.r\");
			msdata=read.table(\"ms${recomb}\")\$V1;
			scrmdata=read.table(\"scrm${recomb}\")\$V1;
			ee=ee_seg(${nsam},${r});
			sdv=sd_recomb(${r},${nsam});
			cat(paste(${nsam},${r},
format(ee,digits=4),format(sdv,digits=4),\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),format(sd(msdata)/sqrt(length(msdata)),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),format(sd(scrmdata)/sqrt(length(scrmdata)),digits=4),
sep=\"\t\"),file=\"${compareRECOMB}\",append=TRUE);cat(\"\n\",file=\"${compareRECOMB}\",append=TRUE);" > dummy.r
			R CMD BATCH dummy.r
		#	rm ms${out} ms${Trees} ms${tmrca} scrm${tmrca} 
			done
		done
	#done

	
	

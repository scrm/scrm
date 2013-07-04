#!/bin/bash

## compare times of genealogy changes
msr=(10 20 10 50)
rep=10
seqlen=100000
msNsample=(2 3 7 10)
compareRECOMB=compareRECOMB
rm ${compareRECOMB}
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
			csplit ms${out} -n 1 '/^///' '{*}'
			rm ms${recomb}
			for file in $(seq 1 1 ${rep})
				do 
				cat xx${file} | grep ";" | wc -l >> ms${recomb}
				done
			scrm ${nsam} ${rep} -r ${r} ${seqlen} -T | tail -n +4  > scrm${out}
			
			rm xx*
			csplit scrm${out} -n 1  '/^///' '{*}'
			rm scrm${recomb}
			for file in $(seq 1 1 ${rep})
				do 
				cat xx${file} | grep ";" | wc -l >> scrm${recomb}
				done
			echo "rm(list=ls());msdata=read.table(\"ms${recomb}\")\$V1;scrmdata=read.table(\"scrm${recomb}\")\$V1;cat(paste(${nsam},${t},mean(msdata),sd(msdata),sd(msdata)/sqrt(length(msdata)),mean(scrmdata),sd(scrmdata),sd(scrmdata)/sqrt(length(scrmdata)),sep=\"\t\"),file=\"${compareRECOMB}\",append=TRUE);cat(\"\n\",file=\"${compareRECOMB}\",append=TRUE);" > dummy.r
			R CMD BATCH dummy.r
		#	rm ms${out} ms${Trees} ms${tmrca} scrm${tmrca} 
			done
		done
	#done

	
	

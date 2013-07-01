#!/bin/bash

#Compare summary statistics of ms and scrm for TMRCA, number of mutation, and number of recombination

msNsample=(2 3 4 7 10 20)

rep=100

## compare TMRCA
compareTMRCA=compareTMRCA
rm ${compareTMRCA}
theta=10
npop=1000000
for nsam in "${msNsample[@]}"
	do
	prefix=${nsam}sample
	out=${prefix}out
	tmrca=${prefix}tmrca
	Trees=${prefix}Trees
	ms ${nsam} ${rep} -t ${theta} -T | tail -n +4 | grep -v "//" > ms${out}
	cat ms${out} | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}
	hybrid-Lambda -gt ms${Trees} -tmrca ms${tmrca}
	scrm ${nsam} ${rep} -t ${theta} -npop ${npop} -tmrca scrm${tmrca}
	echo "rm(list=ls());msdata=read.table(\"ms${tmrca}\")\$V1*4*${npop};scrmdata=read.table(\"scrm${tmrca}\")\$V1;cat(paste(${nsam},mean(msdata),sd(msdata),sd(msdata)/sqrt(length(msdata)),mean(scrmdata),sd(scrmdata),sd(scrmdata)/sqrt(length(scrmdata)),sep=\"\t\"),file=\"${compareTMRCA}\",append=TRUE);cat(\"\n\",file=\"${compareTMRCA}\",append=TRUE);" > dummy.r
	R CMD BATCH dummy.r
	rm ms${out} ms${Trees} ms${tmrca} scrm${tmrca} 
	done


###########################

mst=(10 20 50 100)
msNsample=(2 3 7 10)

npop=20000

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
	
		scrm ${nsam} ${rep} -t ${t} -npop ${npop} -T scrm${out}
		cat scrm${out} | grep "segsites" | sed -e "s/segsites: //" > scrm${nseg}
		
		echo "rm(list=ls());msdata=read.table(\"ms${nseg}\")\$V1;scrmdata=read.table(\"scrm${nseg}\")\$V1;cat(paste(${nsam},${t},mean(msdata),sd(msdata),sd(msdata)/sqrt(length(msdata)),mean(scrmdata),sd(scrmdata),sd(scrmdata)/sqrt(length(scrmdata)),sep=\"\t\"),file=\"${compareSEG}\",append=TRUE);cat(\"\n\",file=\"${compareSEG}\",append=TRUE);" > dummy.r
		R CMD BATCH dummy.r
		rm ms${out} ms${Trees} ms${tmrca}  scrm${out} scrm${tmrca} ms${nseg}  scrm${nseg}
		done
	done


## compare times of genealogy changes
msr=(10 20 10 50)
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
			scrm ${nsam} ${rep} -r ${r} ${seqlen} -npop ${npop} -T scrm${out}
			
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

	
	

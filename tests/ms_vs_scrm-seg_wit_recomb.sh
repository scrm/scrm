#!/bin/bash

## compare number of mutations with recombinations
# this pairwise_test.sh is a sub test of this 

mkdir test-seg-recomb
cd test-seg-recomb
rm *pdf

msr=(11 21 10 50 100)
mst=(10 20 50 100 10)
#msr=(100)
#mst=(25)

rep=100000
seqlen=100000
msNsample=(2 3 7 10)
compareSEG=compareSEG-RECOMB
echo -e "compare number of mutations with recombinations for ${rep} replicates 
\t\t\t|\t\t\t|\t ms \t\t|\t scrm\t\t|\tTest1\t\t\t|\tTest2\t\t\t
Nsam\ttheta\trho\t|\tmean\tstdv\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value\t\t|\tstats\tp-value\t\t" > ${compareSEG}


for t in "${mst[@]}"
	do
	for r in "${msr[@]}"
		do
		for nsam in "${msNsample[@]}"
			do
			prefix=${nsam}sample${r}rho${t}theta
			nseg=${prefix}NumOfSeg
			ms ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} | tail -n +4 | grep -v "//" | grep "segsites" | sed -e "s/segsites: //" > ms${nseg}
	
			scrm ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} | tail -n +4 | grep -v "//"  | grep "segsites" | sed -e "s/segsites: //" > scrm${nseg}
			
			echo "rm(list=ls());
		source(\"../fun_src.r\");
		msdata=read.table(\"ms${nseg}\")\$V1;
		scrmdata=read.table(\"scrm${nseg}\")\$V1;
			ee=ee_seg(${nsam},${t});
			sdv=sd_seg_recomb(${nsam},${t},${r});
			datamax=max(msdata,scrmdata);
		mstable=table(factor(msdata,levels=1:datamax))
		scrmtable=table(factor(scrmdata,levels=1:datamax))
pdf(\"${nsam}sample${t}mut${r}recomb_NUM_MUT.pdf\");
#plot(mstable,col=\"red\",ylab=\"Frequency\",xlab=\"Number of mutations\");
#lines(scrmtable,col=\"blue\");
plot(as.numeric(names(mstable)), mstable/length(msdata),pch=16,col=\"red\",ylab=\"Frequency\",xlab=\"Number of mutations\");
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


cat(paste(${nsam},${t},${r},\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"|\",
format(combined_scrm_ms_test\$statistic,digits=4),format(combined_scrm_ms_test\$p.value,scientific = TRUE),\"|\",
format(combined_ms_scrm_test\$statistic,digits=4),format(combined_ms_scrm_test\$p.value,scientific = TRUE),
sep=\"\t\"),file=\"${compareSEG}\",append=TRUE);cat(\"\n\",file=\"${compareSEG}\",append=TRUE);
" > dummy.r
			R CMD BATCH dummy.r
			rm ms${nseg}  scrm${nseg}
			done
		done
	done

	
	

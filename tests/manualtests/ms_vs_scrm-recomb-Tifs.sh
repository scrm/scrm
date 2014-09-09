#!/bin/bash

## compare times of genealogy changes
mkdir test-recomb-Tifs
cd test-recomb-Tifs
rm *pdf

msr=(10 20 10 50)
msNsample=(2 3 7 10 20)

rep=100000
seqlen=100000

compareRECOMB=compareRECOMB
echo -e "compare number of recombination for ${rep} replicates 
\t\t|\t\t\t|\tms \t\t|\t scrm\t\t|\tTest1\t\t\t|\tTest2\t\t\t
Nsam\trho\t|\tmean\tstdv\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value\t\t|\tstats\tp-value\t\t" >${compareRECOMB}


#for t in "${mst[@]}"
#	do
	for r in "${msr[@]}"
		do
		for nsam in "${msNsample[@]}"
			do
			prefix=${nsam}sample${r}recomb
			#out=${prefix}out
			recomb=${prefix}Recomb
			ms ${nsam} ${rep} -r ${r} ${seqlen} -T | tail -n +4  |	gawk '/^\/\//{f="xx"++d} f{print > f} ' 
			
			for file in $(seq 1 1 ${rep})
				do 
				cat xx${file} | grep ";" | wc -l >> ms${recomb}
				done
				
			find . -name "xx*" -print0 | xargs -0 rm
				
			scrm ${nsam} ${rep} -r ${r} ${seqlen} -Tifs | tail -n +4  | gawk '/^\/\//{f="xx"++d} f{print > f} ' 
			
			for file in $(seq 1 1 ${rep})
				do 
				cat xx${file} | grep ";" | wc -l >> scrm${recomb}
				done
			find . -name "xx*" -print0 | xargs -0 rm	
				
			echo "rm(list=ls());
			source(\"../fun_src.r\");
			msdata=read.table(\"ms${recomb}\")\$V1;
			scrmdata=read.table(\"scrm${recomb}\")\$V1;
			ee=ee_seg(${nsam},${r});
			sdv=sd_recomb(${r},${nsam});
datamax=max(msdata,scrmdata);
		mstable=table(factor(msdata,levels=1:datamax))
		scrmtable=table(factor(scrmdata,levels=1:datamax))
pdf(\"${nsam}sample${r}recomb_NUM_RECOMB.pdf\");
#plot(mstable,col=\"red\",ylab=\"Frequency\",xlab=\"Number of recombinations\");
#lines(scrmtable,col=\"blue\");
plot(as.numeric(names(mstable)), mstable/length(msdata),pch=16,col=\"red\",ylab=\"Frequency\",xlab=\"Number of recombinations\");
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

cat(paste(${nsam},${r},\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"|\",
format(combined_scrm_ms_test\$statistic,digits=4),format(combined_scrm_ms_test\$p.value,scientific = TRUE),\"|\",
format(combined_ms_scrm_test\$statistic,digits=4),format(combined_ms_scrm_test\$p.value,scientific = TRUE),
sep=\"\t\"),file=\"${compareRECOMB}\",append=TRUE);cat(\"\n\",file=\"${compareRECOMB}\",append=TRUE);
" > dummy.r
			R CMD BATCH dummy.r
			rm ms${recomb} scrm${recomb} 
			done
		done
	#done

	
	

#!/bin/bash

###########################

mkdir test-seg
cd test-seg
rm *pdf

mst=(10 20 50 100)
msNsample=(3 7 10)

#mst=(10)
#msNsample=(2)

rep=10000
#npop=20000

## compare number of segregating sites
compareSEG=compareSEG
echo -e "compare number of mutations for ${rep} replicates \n\t\t|\t\t\t|\tms \t\t|\t scrm\t\t|\tTest1\t\t\t|\tTest2\t\t\t
Nsam\ttheta\t|\tmean\tstdv\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value\t\t|\tstats\tp-value\t\t" >${compareSEG}

for t in "${mst[@]}"
	do
	for nsam in "${msNsample[@]}"
		do
		prefix=${nsam}sample${t}mut
		#out=${prefix}out
		nseg=${prefix}NumOfSeg
		ms ${nsam} ${rep} -t ${t} -T | tail -n +4 | grep -v "//" | grep "segsites" | sed -e "s/segsites: //" > ms${nseg}

		scrm ${nsam} ${rep} -t ${t} | tail -n +4 | grep -v "//"  | grep "segsites" | sed -e "s/segsites: //" > scrm${nseg}

		echo "rm(list=ls());
		source(\"../fun_src.r\");
		msdata=read.table(\"ms${nseg}\")\$V1;
		scrmdata=read.table(\"scrm${nseg}\")\$V1;
		ee=ee_seg(${nsam},${t});
		sdv=sd_seg_norecomb(${nsam},${t});
		datamax=max(msdata,scrmdata);
		mstable=table(factor(msdata,levels=1:datamax))
		scrmtable=table(factor(scrmdata,levels=1:datamax))
pdf(\"${nsam}sample${t}mut_NUM_MUT.pdf\");
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

cat(paste(${nsam},${t},\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
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





#chistat=0
#for (i in (1:dim(combinedtable)[1])){
	#for (j in (1:dim(combinedtable)[2])){
		#E=sum(combinedtable[i,])*sum(combinedtable[,j])/sum(combinedtable)
#sum(combinedtable[i,])*sum(combinedtable[,j])/sum(combinedtable)
		#chistat = chistat + (E-combinedtable[i,j])^2/E
	#}

#}

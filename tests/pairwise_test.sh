#!/bin/bash

## compare times of genealogy changes
msr=(0 10 100)

rep=100000
seqlen=100000
#msNsample=(2 3 7 10)
nsam=2
compareRECOMB=compareSEG-RECOMB
#echo -e "compare number of mutations for ${rep} replicates \n\t\t\t|\t\t\t|\t\t ms \t\t|\t\t scrm\nNsam\ttheta\trho\t|\tmean\tstdv\t|\tmean\tstdv\tstd err\t|\tmean\tstdv \tstd err" > ${compareRECOMB}

t=10

echo "rm(list=ls());
pdf(\"pairwise_diff.pdf\");
plot(c(0,50),c(0,0.05),type=\"n\")
for (r in c(0,10,100)){
	ms=read.table(paste(\"ms${nsam}sample\",r,\"rho${t}thetasegRecomb\",sep=\"\"));
	scrm=read.table(paste(\"scrm${nsam}sample\",r,\"rho${t}thetasegRecomb\",sep=\"\"));
	mstable=table(ms)
	scrmtable=table(scrm)
	lines(as.numeric(names(mstable)), mstable/dim(ms)[1],pch=\".\",col=\"red\")
	lines(as.numeric(names(scrmtable)), scrmtable/dim(scrm)[1],pch=\".\",col=\"blue\")
	}
dev.off()
" > dummy.r


for r in "${msr[@]}"
	do
	#for nsam in "${msNsample[@]}"
		#do
		prefix=${nsam}sample${r}rho${t}theta
		out=${prefix}out
		segrecomb=${prefix}segRecomb

		ms ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} -T | tail -n +4 | grep -v "//" > ms${out}
		cat ms${out} | grep "segsites" | sed -e "s/segsites: //" > ms${segrecomb}

		scrm ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} | tail -n +4 | grep -v "//" > scrm${out}
		cat scrm${out} | grep "segsites" | sed -e "s/segsites: //" > scrm${segrecomb}

		
		#echo "rm(list=ls());
		#source(\"fun_src.r\");
		#msdata=read.table(\"ms${segrecomb}\")\$V1;
		#scrmdata=read.table(\"scrm${segrecomb}\")\$V1;
		#ee=ee_seg(${nsam},${r});
		#sdv=sd_seg_recomb(${nsam},${t},${r});
		#cat(paste(${nsam},${t},${r},\"|\",
#format(ee,digits=4),format(sdv,digits=4),\"|\",
#format(mean(msdata),digits=4),format(sd(msdata),digits=4),format(sd(msdata)/sqrt(length(msdata)),digits=4),\"|\",
#format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),format(sd(scrmdata)/sqrt(length(scrmdata)),digits=4),
#sep=\"\t\"),file=\"${compareRECOMB}\",append=TRUE);cat(\"\n\",file=\"${compareRECOMB}\",append=TRUE);" > dummy.r
		#R CMD BATCH dummy.r
		#rm ms${out} ms${segrecomb} scrm${out} scrm${segrecomb} 
		#done
	done
	

R CMD BATCH dummy.r

	
	

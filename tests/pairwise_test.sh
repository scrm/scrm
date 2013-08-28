#!/bin/bash

## Need to review this again... is PSk the probability of observing k number of mutations? (Hein 2005 2.30, Wakeley 2008 4.3) when there is no recombination???
msr=(0 10 100)
msr=(0)
rep=100000
seqlen=100000
nsam=2
compareRECOMB=compareSEG-RECOMB

t=10

echo "rm(list=ls());
pdf(\"pairwise_diff.pdf\");
plot(c(0,50),c(0,0.1),type=\"n\")
for (r in c(0,10,100)){
	ms=read.table(paste(\"ms${nsam}sample\",r,\"rho${t}thetasegRecomb\",sep=\"\"));
	scrm=read.table(paste(\"scrm${nsam}sample\",r,\"rho${t}thetasegRecomb\",sep=\"\"));
	mstable=table(ms)
	scrmtable=table(scrm)
	points(as.numeric(names(mstable)), mstable/dim(ms)[1],pch=\".\",col=\"red\")
	points(as.numeric(names(scrmtable)), scrmtable/dim(scrm)[1],pch=\".\",col=\"blue\")
	}
dev.off()
" > dummy.r


for r in "${msr[@]}"
	do
		prefix=${nsam}sample${r}rho${t}theta
		out=${prefix}out
		segrecomb=${prefix}segRecomb

		ms ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} -T | tail -n +4 | grep -v "//" | grep "segsites" | sed -e "s/segsites: //" > ms${segrecomb}

		scrm ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} | tail -n +4 | grep -v "//"  | grep "segsites" | sed -e "s/segsites: //" > scrm${segrecomb}

	done
	

R CMD BATCH dummy.r

	
	

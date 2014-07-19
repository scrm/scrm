#!/bin/bash
#Compare summary statistics of scrm simulation for the initial genealogy and the final genealogy for TMRCA with recombination events

mkdir test-scrm-first_vs_last
cd test-scrm-first_vs_last
rm *pdf

compareFirstLast=compareFirstLast
rm ${compareFirstLast}


echo "rm(list=ls());
currentcase=scan(\"current_case\",what=\"\");
figuretitle=scan(\"figuretitle\",what=\"\");

firstdata=read.table(\"scrmfirsttmrca\")\$V1;
lastdata=read.table(\"scrmlasttmrca\")\$V1;
test=ks.test(firstdata,lastdata)
pdf(paste(currentcase,figuretitle,\".pdf\",sep=\"\"));
plot(ecdf(firstdata), xlim=range(c(firstdata, lastdata)),col=\"red\", main=currentcase)
plot(ecdf(lastdata), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,digits=4),sep=\"\")))
legend(\"topleft\",c(\"first\",\"last\"), col=c(\"red\",\"blue\"), pch=16)
dev.off();
cat(paste(currentcase,figuretitle,\"|\",
format(mean(firstdata),digits=4),format(sd(firstdata),digits=4),\"|\",
format(mean(lastdata),digits=4),format(sd(lastdata),digits=4),\"|\",test\$statistic,format(test\$p.value,digits=4), 
sep=\"\t\"),file=\"${compareFirstLast}\",append=TRUE);cat(\"\n\",file=\"${compareFirstLast}\",append=TRUE);" > tmrca.r

seqlen=100000

rep=100000

scrmNsample=(2 4 6)
scrmr=(10 20)
scrmt=(10 5)

rm current_case scrm* figuretitle
for t in "${scrmt[@]}"
	do
	for r in "${scrmr[@]}"
		do
		for nsam in "${scrmNsample[@]}"
			do
			rm current_case scrm*Trees scrm*tmrca figuretitle
			prefix=scrm${nsam}sample${r}rho${t}theta
			echo ${prefix} > current_case
			cat current_case

			scrm ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} -T | tail -n +4 | gawk '/^\/\//{f="xx"++d} f{print > f} '
			for file in $(seq 1 1 ${rep})
				do 
				grep ";" xx${file} | sed -e 's/\[.*\]//g' | head -1 >> scrmfirstTrees
				grep ";" xx${file} | sed -e 's/\[.*\]//g' | tail -1 >> scrmlastTrees
				done
			find . -name "xx*" -print0 | xargs -0 rm

			hybrid-Lambda -gt scrmfirstTrees -tmrca scrmfirsttmrca
			hybrid-Lambda -gt scrmlastTrees -tmrca scrmlasttmrca
			echo "TMRCA" > figuretitle
			R CMD BATCH tmrca.r
			rm scrm*tmrca
			hybrid-Lambda -gt scrmfirstTrees -bl scrmfirsttmrca
			hybrid-Lambda -gt scrmlastTrees -bl scrmlasttmrca
			echo "BL" > figuretitle
			R CMD BATCH tmrca.r
			
			done
		done
	done

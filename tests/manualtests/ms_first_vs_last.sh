#!/bin/bash
#Compare summary statistics of ms simulation for the initial genealogy and the final genealogy for TMRCA with recombination events

mkdir test-ms-first_vs_last
cd test-ms-first_vs_last
rm *pdf

compareFirstLast=compareFirstLast
rm ${compareFirstLast}


echo "rm(list=ls());
currentcase=scan(\"current_case\",what=\"\");
figuretitle=scan(\"figuretitle\",what=\"\");

firstdata=read.table(\"msfirsttmrca\")\$V1;
lastdata=read.table(\"mslasttmrca\")\$V1;
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

msNsample=(2 4 6)
msr=(10 20)
mst=(10 5)

for t in "${mst[@]}"
	do
	for r in "${msr[@]}"
		do
		for nsam in "${msNsample[@]}"
			do
			rm current_case ms*Trees ms*tmrca figuretitle
			prefix=ms${nsam}sample${r}rho${t}theta
			echo ${prefix}
			echo "${prefix}" > current_case
			#rm ms* scrm*
			ms ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} -T | tail -n +4 | gawk '/^\/\//{f="xx"++d} f{print > f} '
			for file in $(seq 1 1 ${rep})
				do 
				grep ";" xx${file} | sed -e 's/\[.*\]//g' | head -1 >> msfirstTrees
				grep ";" xx${file} | sed -e 's/\[.*\]//g' | tail -1 >> mslastTrees
				done
				#hybrid-Lambda -gt msfirstTrees -bl msfirstbl
			find . -name "xx*" -print0 | xargs -0 rm
			hybrid-Lambda -gt msfirstTrees -tmrca msfirsttmrca
			hybrid-Lambda -gt mslastTrees -tmrca mslasttmrca
			echo "TMRCA" > figuretitle
			R CMD BATCH tmrca.r
			rm ms*tmrca
			hybrid-Lambda -gt msfirstTrees -bl msfirsttmrca
			hybrid-Lambda -gt mslastTrees -bl mslasttmrca
			echo "BL" > figuretitle
			R CMD BATCH tmrca.r
			
			done
		done
	done

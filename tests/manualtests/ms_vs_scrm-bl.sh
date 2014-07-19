#!/bin/bash

#Compare summary statistics of ms and scrm for BL, number of mutation, and number of recombination

mkdir test-bl
cd test-bl
rm *pdf

msNsample=(2 3 4 7 10 20)

rep=100000

## compare BL
compareBL=compareBL
rm ${compareBL}
echo -e "compare BL for ${rep} replicates \n\t|Theoretical\t\t|\t ms \t\t|\t scrm\t\t|\tKS test\nNsam\t|\tmean\tstdv\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value" >${compareBL}
#theta=10

npop=1000000
for nsam in "${msNsample[@]}"
	do
	prefix=${nsam}sample
	out=${prefix}out
	bl=${prefix}bl
	Trees=${prefix}Trees
	ms ${nsam} ${rep} -T | tail -n +4 | grep -v "//" > ms${out}
	cat ms${out} | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}
	hybrid-Lambda -gt ms${Trees} -bl ms${bl}
	scrm ${nsam} ${rep} -T | tail -n +4 | grep -v "//" > scrm${out}
	cat scrm${out} | grep ";" | sed -e 's/\[.*\]//g' > scrm${Trees}
	hybrid-Lambda -gt scrm${Trees} -bl scrm${bl}
	echo "rm(list=ls());
source(\"../fun_src.r\");
msdata=read.table(\"ms${bl}\")\$V1;
scrmdata=read.table(\"scrm${bl}\")\$V1;
ee=ee_bl(${nsam});
sdv=sd_bl(${nsam});
test=ks.test(msdata,scrmdata)
pdf(paste($nsam,\"sampleBL-KStest.pdf\",sep=\"\"));
plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=paste($nsam,\"sample BL KS test\",sep=\"\"))
plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,digits=4),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
dev.off();
cat(paste(${nsam},\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"|\",test\$statistic,format(test\$p.value,digits=4), 
sep=\"\t\"),file=\"${compareBL}\",append=TRUE);cat(\"\n\",file=\"${compareBL}\",append=TRUE);" > dummy.r
	R CMD BATCH dummy.r
	rm ms${out} ms${Trees} ms${bl} scrm${out} scrm${Trees} scrm${bl} 
	done


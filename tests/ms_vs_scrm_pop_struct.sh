#!/bin/bash

mkdir test-POP
cd test-POP
rm *pdf


rep=100000

## compare population sturture for a single population data
comparePop=comparePop
rm ${comparePop}

theta=10


echo "rm(list=ls());
#source(\"../fun_src.r\");
figuretitle=scan(\"figuretitle\",what=\"\");
currentcase=scan(\"current_case\",what=\"\");
msdata=read.table(paste(\"ms\",\"data\",sep=\"\"))\$V1;
scrmdata=read.table(paste(\"scrm\",\"data\",sep=\"\"))\$V1;
#ee=1#ee_tmrca(${nsam});
#sdv=1#sd_tmrca(${nsam});
test=ks.test(msdata,scrmdata)
pdf(paste(currentcase,figuretitle,\".pdf\",sep=\"\"));
plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=currentcase)
plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,digits=4),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
dev.off();
#cat(paste(currentcase,\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
cat(paste(currentcase, figuretitle , \"\n\",\"|\",
format(mean(msdata),scientific = TRUE),format(sd(msdata),scientific = TRUE),\"|\",
format(mean(scrmdata),scientific = TRUE),format(sd(scrmdata),scientific = TRUE),\"|\",format(test\$statistic,scientific = TRUE),format(test\$p.value,scientific = TRUE), 
sep=\"\t\"),file=\"${comparePop}\",append=TRUE);cat(\"\n\",file=\"${comparePop}\",append=TRUE);" > dummy.r

foo(){
	cut -f 2 ms_stats > msdata
	cut -f 2 scrm_stats > scrmdata
	echo "Pairewise_difference" > figuretitle
	R CMD BATCH dummy.r
	
	cut -f 6 ms_stats > msdata
	cut -f 6 scrm_stats > scrmdata
	echo "Tajima_D" > figuretitle
	R CMD BATCH dummy.r
	
	cut -f 8 ms_stats > msdata
	cut -f 8 scrm_stats > scrmdata
	echo "theta_H" > figuretitle
	R CMD BATCH dummy.r
	
	cut -f 10 ms_stats > msdata
	cut -f 10 scrm_stats > scrmdata
	echo "H" > figuretitle
	R CMD BATCH dummy.r
	}

#case 1 
echo "4_samples" > current_case

ms 4 ${rep} -t ${theta} -eN 0.8 0.1 -eN 2.0 4.0 | sample_stats > ms_stats
scrm 4 ${rep} -t ${theta} -eN 0.8 0.1 -eN 2.0 4.0  | sample_stats > scrm_stats

foo


#case 2

echo "5_samples" > current_case

ms 5 ${rep} -t ${theta} -eN 0.8 10.1 -eN 2.0 10.0 | sample_stats > ms_stats
scrm 5 ${rep} -t ${theta} -eN 0.8 10.1 -eN 2.0 10.0 | sample_stats > scrm_stats

foo

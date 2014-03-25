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
		#ee=ee_seg(${nsam},${t});
		#sdv=sd_seg_norecomb(${nsam},${t});
		datamax=max(msdata,scrmdata);
		mstable=table(factor(msdata,levels=1:datamax))
		scrmtable=table(factor(scrmdata,levels=1:datamax))
pdf(paste(currentcase,figuretitle,\".pdf\",sep=\"\"));
plot(as.numeric(names(mstable)), mstable/length(msdata),pch=16,col=\"red\",ylab=\"Frequency\",xlab=figuretitle);
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

#cat(paste(${nsam},${t},\"|\",format(ee,digits=4),format(sdv,digits=4),\"|\",
cat(paste(currentcase, figuretitle , \"\n\",\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"||\",
format(combined_scrm_ms_test\$statistic,digits=4),format(combined_scrm_ms_test\$p.value,scientific = TRUE),\"||\",
format(combined_ms_scrm_test\$statistic,digits=4),format(combined_ms_scrm_test\$p.value,scientific = TRUE),
sep=\"\t\"),file=\"${comparePop}\",append=TRUE);cat(\"\n\",file=\"${comparePop}\",append=TRUE);
" > chisq.r

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
format(mean(msdata),scientific = TRUE),format(sd(msdata),scientific = TRUE),\"||\",
format(mean(scrmdata),scientific = TRUE),format(sd(scrmdata),scientific = TRUE),\"|\",format(test\$statistic,scientific = TRUE),format(test\$p.value,scientific = TRUE), 
sep=\"\t\"),file=\"${comparePop}\",append=TRUE);cat(\"\n\",file=\"${comparePop}\",append=TRUE);" > ks.r

echo "rm(list=ls());
#source(\"../fun_src.r\");
figuretitle=scan(\"figuretitle\",what=\"\");
currentcase=scan(\"current_case\",what=\"\");
msdata=read.table(\"mstmrca\")\$V1;
scrmdata=read.table(\"scrmtmrca\")\$V1;
#ee=ee_tmrca(${nsam});
#sdv=sd_tmrca(${nsam});
test=ks.test(msdata,scrmdata)
pdf(paste(currentcase,figuretitle,\".pdf\",sep=\"\"));
plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=currentcase)
plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,digits=4),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
dev.off();
cat(paste(currentcase,figuretitle , \"\n\",\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"|\",test\$statistic,format(test\$p.value,digits=4), 
sep=\"\t\"),file=\"${comparePop}\",append=TRUE);cat(\"\n\",file=\"${comparePop}\",append=TRUE);" > tmrca.r

echo "rm(list=ls());
#source(\"../fun_src.r\");
figuretitle=scan(\"figuretitle\",what=\"\");
currentcase=scan(\"current_case\",what=\"\");
msdata=read.table(\"msbl\")\$V1;
scrmdata=read.table(\"scrmbl\")\$V1;
test=ks.test(msdata,scrmdata)
pdf(paste(currentcase,figuretitle,\".pdf\",sep=\"\"));
plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=currentcase)
plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,digits=4),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
dev.off();
cat(paste(currentcase,figuretitle , \"\n\",\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"|\",test\$statistic,format(test\$p.value,digits=4), 
sep=\"\t\"),file=\"${comparePop}\",append=TRUE);cat(\"\n\",file=\"${comparePop}\",append=TRUE);" > bl.r

#format(ee,digits=4),format(sdv,digits=4),\"|\",
foo(){
	echo "TMRCA" > figuretitle
	R CMD BATCH tmrca.r

	echo "BL" > figuretitle
	R CMD BATCH bl.r

	cut -f 6 ms_stats > msdata
	cut -f 6 scrm_stats > scrmdata
	echo "Tajima_D" > figuretitle
	R CMD BATCH ks.r

	cut -f 2 ms_stats > msdata
	cut -f 2 scrm_stats > scrmdata
	echo "Pairewise_difference" > figuretitle
	R CMD BATCH chisq.r

	cut -f 8 ms_stats > msdata
	cut -f 8 scrm_stats > scrmdata
	echo "theta_H" > figuretitle
	R CMD BATCH chisq.r

	cut -f 10 ms_stats > msdata
	cut -f 10 scrm_stats > scrmdata
	echo "H" > figuretitle
	R CMD BATCH chisq.r
	}

#case 1 
echo "2_samples_case1" > current_case
rm ms* scrm*
ms 2 ${rep} -t ${theta} -eN 0.4 10.01 -eN 1 0.01 -T > msout
scrm 2 ${rep} -t ${theta} -eN 0.4 10.01 -eN 1 0.01  -T > scrmout


cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 1 
echo "2_samples_case1.1" > current_case
rm ms* scrm*
ms 2 ${rep} -t ${theta} -eN 0 1 -eN 0.4 10.01 -eN 1 0.01 -T > msout
scrm 2 ${rep} -t ${theta} -eN 0 1 -eN 0.4 10.01 -eN 1 0.01  -T > scrmout


cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 2
echo "2_samples_case1.2" > current_case
rm ms* scrm*
ms 2 ${rep} -t ${theta} -eN 0 10 -eN 0.4 10.01 -eN 1 0.01 -T > msout
scrm 2 ${rep} -t ${theta} -eN 0 10 -eN 0.4 10.01 -eN 1 0.01  -T > scrmout


cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 2
echo "2_samples_case2" > current_case
rm ms* scrm*
ms 2 ${rep} -t ${theta} -eN .5 0.01 -T > msout
scrm 2 ${rep} -t ${theta} -eN .5 0.01  -T > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo


##case 3

echo "5_samples" > current_case
rm ms* scrm*

ms 5 ${rep} -t ${theta} -eN 0.5 10.0 -T > msout
scrm 5 ${rep} -t ${theta} -eN 0.5 10.0 -T > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 4
echo "4_samples" > current_case
rm ms* scrm*

ms 4 ${rep} -t ${theta} -eN 0.8 15  -T > msout
scrm 4 ${rep} -t ${theta} -eN 0.8 15  -T > scrmout

cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 5
echo "6_samples_case1" > current_case
rm ms* scrm*
ms 6 ${rep} -t ${theta} -eN 1 .1 -eN 3 10 -T > msout
scrm 6 ${rep} -t ${theta} -eN 1 .1 -eN 3 10  -T > scrmout
cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

#case 6
echo "6_samples_case2" > current_case
rm ms* scrm*
ms 6 ${rep} -t ${theta}  -eN .3 10 -T > msout
scrm 6 ${rep} -t ${theta}  -eN .3 10  -T > scrmout
cat msout | sample_stats > ms_stats
cat msout | grep ";" | sed -e 's/\[.*\]//g' > msTrees
hybrid-Lambda -gt msTrees -tmrca mstmrca
hybrid-Lambda -gt msTrees -bl msbl


cat scrmout | sample_stats > scrm_stats
cat scrmout | grep ";" | sed -e 's/\[.*\]//g' > scrmTrees
hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
hybrid-Lambda -gt scrmTrees -bl scrmbl

foo

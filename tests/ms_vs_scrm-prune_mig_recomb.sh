#!/bin/bash

#Compare summary statistics of ms and scrm for migration

mkdir test-prune-mig-pair
cd test-prune-mig-pair
rm *pdf


rep=10000
seqlen=1000000
#msr=(10 20 10 50)
r=100

## compare TMRCA
compareMIG=compareMIG
rm ${comparePop}


#rm ${compareTMRCA}
#rm ${compareBL}
#echo -e "compare TMRCA for ${rep} replicates 
#\t|\t ms \t\t|\t scrm\t\t|\tKS test
#Case\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value" >${compareTMRCA}

#echo -e "compare BL for ${rep} replicates 
#\t|\t ms \t\t|\t scrm\t\t|\tKS test
#Case\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value" >${compareBL}

theta=100

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
sep=\"\t\"),file=\"${compareMIG}\",append=TRUE);cat(\"\n\",file=\"${compareMIG}\",append=TRUE);
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
sep=\"\t\"),file=\"${compareMIG}\",append=TRUE);cat(\"\n\",file=\"${compareMIG}\",append=TRUE);" > ks.r

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
sep=\"\t\"),file=\"${compareMIG}\",append=TRUE);cat(\"\n\",file=\"${compareMIG}\",append=TRUE);" > tmrca.r

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
sep=\"\t\"),file=\"${compareMIG}\",append=TRUE);cat(\"\n\",file=\"${compareMIG}\",append=TRUE);" > bl.r


foo(){
    #cut -f 2 mstime > mstmrca
	#cut -f 2 scrmtime > scrmtmrca
	echo "TMRCA" > figuretitle
	R CMD BATCH tmrca.r

    #cut -f 3 mstime > msbl
	#cut -f 3 scrmtime > scrmbl
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
	
mstime(){
    cat msout | gawk '/^\/\//{f="xx"++d} f{print > f} '
    for file in $(seq 1 1 ${rep})
        do 
        grep ";" xx${file} | sed -e 's/\[.*\]//g' | tail -1 >> msTrees
        done
        hybrid-Lambda -gt msTrees -tmrca mstmrca
        hybrid-Lambda -gt msTrees -bl msbl
        find . -name "xx*" -print0 | xargs -0 rm    
    }	
    
scrmtime(){
    cat scrmout | gawk '/^\/\//{f="xx"++d} f{print > f} '
    for file in $(seq 1 1 ${rep})
        do 
        grep ";" xx${file} | sed -e 's/\[.*\]//g' | tail -1 >> scrmTrees
        done
        hybrid-Lambda -gt scrmTrees -tmrca scrmtmrca
        hybrid-Lambda -gt scrmTrees -bl scrmbl
        find . -name "xx*" -print0 | xargs -0 rm    
    }	
	
	
#case 1 
#2 sub population, 1 sample from each subpopulation, mutation rate is 5
echo "2groups1sam1sam_case1_mig5" > current_case
rm ms* scrm*

#ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}

ms 2 ${rep} -t ${theta} -r ${r} ${seqlen} -I 2 1 1 5.0 -T > msout
mstime 

scrm 2 ${rep} -t ${theta} -r ${r} ${seqlen}  -I 2 1 1 -ma x 5.0 5.0 x -T -l 10000 > scrmout
scrmtime

cat msout | sample_stats > ms_stats

cat scrmout | sample_stats > scrm_stats
#cat scrmout | grep "time:" >  scrmtime

foo

#case 2
#2 sub population, 1 sample from each subpopulation, mutation rate is 5
echo "2groups1sam1sam_case2_mig5" > current_case
rm ms* scrm*

#ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}

ms 2 ${rep} -t ${theta} -r ${r} ${seqlen}  -I 2 1 1 5.0 -T > msout
mstime

scrm 2 ${rep} -t ${theta} -r ${r} ${seqlen}  -I 2 1 1 5.0 -T -L -l 10000> scrmout
scrmtime


cat msout | sample_stats > ms_stats

cat scrmout | sample_stats > scrm_stats
#cat scrmout | grep "time:" >  scrmtime

foo


#case 3
#2 sub population, 1 sample from each subpopulation, mutation rate is 5
echo "2groups1sam1sam_case3_mig5" > current_case
rm ms* scrm*

#ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}

ms 2 ${rep} -t ${theta} -r ${r} ${seqlen}  -I 2 1 1 -ma x 5 0 x -ej 0.5 2 1 -T > msout
mstime

scrm 2 ${rep} -t ${theta} -r ${r} ${seqlen}  -I 2 1 1 -ma x 5 0 x -ej 0.5 2 1 -T -L -l 10000 > scrmout
scrmtime

cat msout | sample_stats > ms_stats

cat scrmout | sample_stats > scrm_stats
#cat scrmout | grep "time:" >  scrmtime

foo

#case 3
#2 sub population, 1 sample from each subpopulation, mutation rate is 5
echo "2groups1sam1sam_case4_mig5" > current_case
rm ms* scrm*

#ms 4 ${rep} -t ${theta} -I 2 2 2 -ma x 5.0 5.0 x -T | tail -n +4 | grep -v "//" | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}

ms 2 ${rep} -t ${theta} -r ${r} ${seqlen}  -I 2 1 1 -ma x 0 5 x -ej 0.5 2 1 -T > msout
mstime

scrm 2 ${rep} -t ${theta} -r ${r} ${seqlen}  -I 2 1 1 -ma x 0 5 x -ej 0.5 2 1 -T -L -l 10000 > scrmout
scrmtime

cat msout | sample_stats > ms_stats

cat scrmout | sample_stats > scrm_stats
#cat scrmout | grep "time:" >  scrmtime

foo


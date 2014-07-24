#!/bin/bash

#Compare summary statistics of ms and scrm for TMRCA with recombination events
mkdir test-tmrca-last
cd test-tmrca-last


seqlen=100000

rep=10000


#echo -e "rm(list=ls());
#for (i in (1:${rep})){
	#p=read.table(paste(\"xx\",i,\"Trees_freq\",sep=\"\"))\$V1/${seqlen}; 
	#T=read.table(paste(\"xx\",i,\"Trees_tmrca\",sep=\"\")); 
	#bl=read.table(paste(\"xx\",i,\"Trees_bl\",sep=\"\")); 
	#cat(paste(sum(T*p),sum(T^2*p),sum(T^3*p),sum(T^4*p),sep=\"\t\") ,file=\"tmrca_moments\",append=TRUE);
	#cat(\"\n\",file=\"tmrca_moments\",append=TRUE);
	#cat(paste(sum(bl*p),sum(bl^2*p),sum(bl^3*p),sum(bl^4*p),sep=\"\t\") ,file=\"bl_moments\",append=TRUE);
	#cat(\"\n\",file=\"bl_moments\",append=TRUE);
	#}" > compute_moments.r

#msNsample=(2 3 4 7 10 20)
#msNsample=(2 4)
#msr=(10 20 50 100)
#mst=(10 20 50 100)

msNsample=(6)
msr=(10)
mst=(10)


## compare TMRCA
compareTMRCA=compareTMRCA-LAST
rm ${compareTMRCA}

compareBL=compareBL-LAST
rm ${compareBL}

rm *pdf
#theta=10
#echo -e "compare TMRCA for ${rep} replicates \n\t\t\t|\t1st\t\t\t|\t2nd\t\t\t|\t3rd\t\t\t|\t4th\t\t\t
#\t\t\t|\tstats\tp-value\t\t|\tstats\tp-value\t\t|\tstats\tp-value\t\t|\tstats\tp-value\t\t" > ${compareTMRCA}

#echo -e "compare total branch length for ${rep} replicates \n\t\t\t|\t1st\t\t\t|\t2nd\t\t\t|\t3rd\t\t\t|\t4th\t\t\t
#\t\t\t|\tstats\tp-value\t\t|\tstats\tp-value\t\t|\tstats\tp-value\t\t|\tstats\tp-value\t\t" > ${compareBL}

echo -e "compare TMRCA for ${rep} replicates 
\t\t\t|\t ms \t\t|\t scrm\t\t|\tKS test
Nsam\ttheta\trho\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value" >${compareTMRCA}

echo -e "compare BL for ${rep} replicates 
\t\t\t|\t ms \t\t|\t scrm\t\t|\tKS test
Nsam\ttheta\trho\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value" >${compareBL}



#rm tmrca_moments
#rm bl_moments


for t in "${mst[@]}"
	do
	for r in "${msr[@]}"
		do
		for nsam in "${msNsample[@]}"
			do
			#rm ms_tmrca_moments
			#rm scrm_tmrca_moments
			
			#rm ms_bl_moments
			#rm scrm_bl_moments
			
			prefix=${nsam}sample${r}rho${t}theta
			out=${prefix}out
			bl=${prefix}bl
			tmrca=${prefix}tmrca
			#segrecomb=${prefix}segRecomb
			rm msTrees
			ms ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} -T | tail -n +4 | gawk '/^\/\//{f="xx"++d} f{print > f} '
			for file in $(seq 1 1 ${rep})
				do 
				grep ";" xx${file} | sed -e 's/\[.*\]//g' | tail -1 >> msTrees
				done
				hybrid-Lambda -gt msTrees -tmrca ms${tmrca}
				hybrid-Lambda -gt msTrees -bl ms${bl}
				find . -name "xx*" -print0 | xargs -0 rm
				
			rm 	scrmTrees
			scrm ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} -T -l 100 | tail -n +4 | gawk '/^\/\//{f="xx"++d} f{print > f} '
			for file in $(seq 1 1 ${rep})
				do 
				grep ";" xx${file} | sed -e 's/\[.*\]//g' | tail -1 >> scrmTrees
				done
				hybrid-Lambda -gt scrmTrees -tmrca scrm${tmrca}
				hybrid-Lambda -gt scrmTrees -bl scrm${bl}
				find . -name "xx*" -print0 | xargs -0 rm

				echo "rm(list=ls());
source(\"../fun_src.r\");
msdata=read.table(\"ms${tmrca}\")\$V1;
scrmdata=read.table(\"scrm${tmrca}\")\$V1;
#ee=ee_tmrca(${nsam});
#sdv=sd_tmrca(${nsam});
test=ks.test(msdata,scrmdata)
pdf(paste($nsam,\"sample${t}mut${r}TMRCA-KStest.pdf\",sep=\"\"));
plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=paste($nsam,\"sample TMRCA KS test\",sep=\"\"))
plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,scientific = TRUE),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
dev.off();
cat(paste(${nsam},${t},${r},\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"|\",format(test\$statistic,digits=4),format(test\$p.value,scientific = TRUE), 
sep=\"\t\"),file=\"${compareTMRCA}\",append=TRUE);cat(\"\n\",file=\"${compareTMRCA}\",append=TRUE);
rm(list=ls());
source(\"../fun_src.r\");
msdata=read.table(\"ms${bl}\")\$V1;
scrmdata=read.table(\"scrm${bl}\")\$V1;
ee=ee_bl(${nsam});
sdv=sd_bl(${nsam});
test=ks.test(msdata,scrmdata)
pdf(paste($nsam,\"sample${t}mut${r}BL-KStest.pdf\",sep=\"\"));
plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col=\"red\", main=paste($nsam,\"sample BL KS test\",sep=\"\"))
plot(ecdf(scrmdata), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,scientific = TRUE),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
dev.off();
cat(paste(${nsam},${t},${r},\"|\",
format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"|\",format(test\$statistic,digits=4),format(test\$p.value,scientific = TRUE), 
sep=\"\t\"),file=\"${compareBL}\",append=TRUE);cat(\"\n\",file=\"${compareBL}\",append=TRUE);" > dummy.r
			R CMD BATCH dummy.r
			#rm scrm${tmrca} scrm${bl} ms${tmrca} ms${bl}
			done
		done
	done
	
	

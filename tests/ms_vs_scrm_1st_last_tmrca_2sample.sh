#!/bin/bash

#Compare summary statistics of ms and scrm for TMRCA with recombination events
mkdir test-first-last
cd test-first-last


seqlen=100000

rep=100000


nsam=2 
#msr=(10 20 50 100)
#mst=(10 20 50 100)

#msNsample=(3)
msr=(10)
mst=(10)


## compare TMRCA
compareTMRCA=compareTMRCA-LAST
rm ${compareTMRCA}

#compareBL=compareBL-LAST
#rm ${compareBL}

#rm *pdf

echo -e "compare TMRCA for ${rep} replicates 
\t\t\t|\t ms \t\t|\t scrm\t\t|\tKS test
Nsam\ttheta\trho\t|\tmean\tstdv\t|\tmean\tstdv \t|\tstats\tp-value" >${compareTMRCA}


#rm tmrca_moments
#rm bl_moments


for t in "${mst[@]}"
	do
	for r in "${msr[@]}"
		do
			prefix=${nsam}sample${r}rho${t}theta
			echo ${prefix}
			out=${prefix}out
			#bl=${prefix}bl
			tmrca=${prefix}tmrca
			
			rm first_ms${tmrca}_file
			rm last_ms${tmrca}_file
		
			ms ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} -T | tail -n +4 | gawk '/^\/\//{f="xx"++d} f{print > f} '
			for file in $(seq 1 1 ${rep})
				do 
				#grep ";" xx${file} | sed -e 's/\[.*\]//g' | tail -1 >> msTrees
				grep ';' xx${file} | sed -e 's/\[.*://g' -e 's/);//g' | head -1 >> first_ms${tmrca}_file
				grep ';' xx${file} | sed -e 's/\[.*://g' -e 's/);//g' | tail -1 >> last_ms${tmrca}_file
				done
				find . -name "xx*" -print0 | xargs -0 rm

			rm first_scrm${tmrca}_file
			rm last_scrm${tmrca}_file
		
			scrm ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} -T -l 100 | tail -n +4 | gawk '/^\/\//{f="xx"++d} f{print > f} '
			for file in $(seq 1 1 ${rep})
				do 
				#grep ";" xx${file} | sed -e 's/\[.*\]//g' | tail -1 >> scrmTrees
				grep ';' xx${file} | sed -e 's/\[.*://g' -e 's/);//g' | head -1 >> first_scrm${tmrca}_file
				grep ';' xx${file} | sed -e 's/\[.*://g' -e 's/);//g' | tail -1 >> last_scrm${tmrca}_file
				done
				find . -name "xx*" -print0 | xargs -0 rm
				

				echo "rm(list=ls());
source(\"../fun_src.r\");
msfirst = read.table(\"first_ms${tmrca}_file\")\$V1;
mslast = read.table(\"last_ms${tmrca}_file\")\$V1;
scrmfirst = read.table(\"first_scrm${tmrca}_file\")\$V1;
scrmlast = read.table(\"last_scrm${tmrca}_file\")\$V1;
#ee=ee_tmrca(${nsam});
#sdv=sd_tmrca(${nsam});
test=ks.test(mslast,scrmlast)
pdf(paste($nsam,\"sample${t}mut${r}TMRCA-KStest.pdf\",sep=\"\"));
plot(ecdf(msfirst), xlim=range(c(msfirst, mslast,scrmfirst,scrmlast)),col=\"red\", main=paste($nsam,\"sample TMRCA KS test\",sep=\"\"))
plot(ecdf(mslast), add=TRUE, lty=\"dashed\", col=\"blue\")
plot(ecdf(scrmfirst), add=TRUE, lty=\"dashed\", col=\"green\")
plot(ecdf(scrmlast), add=TRUE, lty=\"dashed\", col=\"black\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,scientific = TRUE),sep=\"\")))
legend(\"topleft\",c(\"msfirst\",\"mslast\",\"scrmfirst\",\"scrmlast\"), col=c(\"red\",\"blue\",\"green\",\"black\"), pch=16)
dev.off();
pdf(paste($nsam,\"sample${t}mut${r}TMRCA-KStest_zoomin.pdf\",sep=\"\"));
plot(ecdf(msfirst), xlim=c(1,2.5),ylim=c(0.8,1),col=\"red\", main=paste($nsam,\"sample TMRCA KS test\",sep=\"\"))
plot(ecdf(mslast), add=TRUE, lty=\"dashed\", col=\"blue\")
plot(ecdf(scrmfirst), add=TRUE, lty=\"dashed\", col=\"green\")
plot(ecdf(scrmlast), add=TRUE, lty=\"dashed\", col=\"black\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",test\$statistic,sep=\"\"), paste(\"p-value = \",format(test\$p.value,scientific = TRUE),sep=\"\")))
legend(\"topleft\",c(\"msfirst\",\"mslast\",\"scrmfirst\",\"scrmlast\"), col=c(\"red\",\"blue\",\"green\",\"black\"), pch=16)
dev.off();
#cat(paste(${nsam},${t},${r},\"|\",
#format(mean(msdata),digits=4),format(sd(msdata),digits=4),\"|\",
#format(mean(scrmdata),digits=4),format(sd(scrmdata),digits=4),\"|\",format(test\$statistic,digits=4),format(test\$p.value,scientific = TRUE), 
#sep=\"\t\"),file=\"${compareTMRCA}\",append=TRUE);cat(\"\n\",file=\"${compareTMRCA}\",append=TRUE); " > dummy.r
			R CMD BATCH dummy.r
			#rm scrm${tmrca} scrm${bl} ms${tmrca} ms${bl}
		done
	done
	
	

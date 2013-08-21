#!/bin/bash

#Compare summary statistics of ms and scrm for TMRCA with recombination events
mkdir test-tmrca-recomb
cd test-tmrca-recomb


seqlen=100000

rep=10000


echo -e "rm(list=ls());
for (i in (1:${rep})){
	p=read.table(paste(\"xx\",i,\"Trees_freq\",sep=\"\"))\$V1/${seqlen}; 
	T=read.table(paste(\"xx\",i,\"Trees_tmrca\",sep=\"\")); 
	bl=read.table(paste(\"xx\",i,\"Trees_bl\",sep=\"\")); 
	cat(paste(sum(T*p),sum(T^2*p),sum(T^3*p),sum(T^4*p),sep=\"\t\") ,file=\"tmrca_moments\",append=TRUE);
	cat(\"\n\",file=\"tmrca_moments\",append=TRUE);
	cat(paste(sum(bl*p),sum(bl^2*p),sum(bl^3*p),sum(bl^4*p),sep=\"\t\") ,file=\"bl_moments\",append=TRUE);
	cat(\"\n\",file=\"bl_moments\",append=TRUE);
	}" > compute_moments.r

msNsample=(2 3 4 7 10 20)
msr=(10 20 10 50 100)
mst=(10 20 50 100 10)

#msNsample=(3)
#msr=(10)
#mst=(10)


## compare TMRCA
compareTMRCA=compareTMRCA-RECOMB
rm ${compareTMRCA}

compareBL=compareBL-RECOMB
rm ${compareBL}

rm *pdf
#theta=10
echo -e "compare TMRCA for ${rep} replicates \n\t\t\t|\t1st\t\t\t|\t2nd\t\t\t|\t3rd\t\t\t|\t4th\t\t\t
\t\t\t|\tstats\tp-value\t\t|\tstats\tp-value\t\t|\tstats\tp-value\t\t|\tstats\tp-value\t\t" > ${compareTMRCA}

echo -e "compare total branch length for ${rep} replicates \n\t\t\t|\t1st\t\t\t|\t2nd\t\t\t|\t3rd\t\t\t|\t4th\t\t\t
\t\t\t|\tstats\tp-value\t\t|\tstats\tp-value\t\t|\tstats\tp-value\t\t|\tstats\tp-value\t\t" > ${compareBL}
rm tmrca_moments
rm bl_moments


for t in "${mst[@]}"
	do
	for r in "${msr[@]}"
		do
		for nsam in "${msNsample[@]}"
			do
			rm ms_tmrca_moments
			rm scrm_tmrca_moments
			
			rm ms_bl_moments
			rm scrm_bl_moments
			
			prefix=${nsam}sample${r}rho${t}theta
			out=${prefix}out
			segrecomb=${prefix}segRecomb

			ms ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} -T | tail -n +4 | gawk '/^\/\//{f="xx"++d} f{print > f} '
			for file in $(seq 1 1 ${rep})
				do 
				grep ";" xx${file} | sed -e 's/\[.*\]//g' > xxTrees
				hybrid-Lambda -gt xxTrees -tmrca xx${file}Trees_tmrca
				hybrid-Lambda -gt xxTrees -bl xx${file}Trees_bl
				grep ";" xx${file} | sed -e 's/\[//g' -e 's/\].*;//g' > xx${file}Trees_freq

				done
				R CMD BATCH compute_moments.r
				find . -name "xx*" -print0 | xargs -0 rm
				mv tmrca_moments ms_tmrca_moments
				mv bl_moments ms_bl_moments
				
				
			scrm ${nsam} ${rep} -t ${t} -r ${r} ${seqlen} -T | tail -n +4 | gawk '/^\/\//{f="xx"++d} f{print > f} '
			for file in $(seq 1 1 ${rep})
				do 
				grep ";" xx${file} | sed -e 's/\[.*\]//g' > xxTrees
				hybrid-Lambda -gt xxTrees -tmrca xx${file}Trees_tmrca
				hybrid-Lambda -gt xxTrees -bl xx${file}Trees_bl
				grep ";" xx${file} | sed -e 's/\[//g' -e 's/\].*;//g' > xx${file}Trees_freq

				done
				R CMD BATCH compute_moments.r
				find . -name "xx*" -print0 | xargs -0 rm
				mv tmrca_moments scrm_tmrca_moments
				mv bl_moments scrm_bl_moments

				echo "rm(list=ls());
msdata=read.table(\"ms_tmrca_moments\");
scrmdata=read.table(\"scrm_tmrca_moments\");
m1test=ks.test(msdata\$V1,scrmdata\$V1)
m2test=ks.test(msdata\$V2,scrmdata\$V2)
m3test=ks.test(msdata\$V3,scrmdata\$V3)
m4test=ks.test(msdata\$V4,scrmdata\$V4)
pdf(paste($nsam,\"sampleTMRCAmut\",${t},\"recomb\",${r},\"-KStest.pdf\",sep=\"\"));
par(mfrow=c(2,2))
plot(ecdf(msdata\$V1), xlim=range(c(msdata\$V1, scrmdata\$V1)),col=\"red\", main=\"1st Moment\")
plot(ecdf(scrmdata\$V1), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",m1test\$statistic,sep=\"\"), paste(\"p-value = \",format(m1test\$p.value,scientific = TRUE),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
####
plot(ecdf(msdata\$V2), xlim=range(c(msdata\$V2, scrmdata\$V2)),col=\"red\", main=\"2nd Moment\")
plot(ecdf(scrmdata\$V2), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",m2test\$statistic,sep=\"\"), paste(\"p-value = \",format(m2test\$p.value,scientific = TRUE),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
####
plot(ecdf(msdata\$V3), xlim=range(c(msdata\$V3, scrmdata\$V3)),col=\"red\", main=\"3rd Moment\")
plot(ecdf(scrmdata\$V3), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",m3test\$statistic,sep=\"\"), paste(\"p-value = \",format(m3test\$p.value,scientific = TRUE),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
####
plot(ecdf(msdata\$V4), xlim=range(c(msdata\$V4, scrmdata\$V4)),col=\"red\", main=\"4th Moment\")
plot(ecdf(scrmdata\$V4), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",m4test\$statistic,sep=\"\"), paste(\"p-value = \",format(m4test\$p.value,scientific = TRUE),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
####
dev.off();
cat(paste(${nsam},${t},${r},\"|\",
format(m1test\$statistic,digits=4),format(m1test\$p.value,scientific = TRUE),\"|\",
format(m2test\$statistic,digits=4),format(m2test\$p.value,scientific = TRUE),\"|\",
format(m3test\$statistic,digits=4),format(m3test\$p.value,scientific = TRUE),\"|\",
format(m4test\$statistic,digits=4),format(m4test\$p.value,scientific = TRUE),
sep=\"\t\"),file=\"${compareTMRCA}\",append=TRUE);cat(\"\n\",file=\"${compareTMRCA}\",append=TRUE);

rm(list=ls());
msdata=read.table(\"ms_bl_moments\");
scrmdata=read.table(\"scrm_bl_moments\");
m1test=ks.test(msdata\$V1,scrmdata\$V1)
m2test=ks.test(msdata\$V2,scrmdata\$V2)
m3test=ks.test(msdata\$V3,scrmdata\$V3)
m4test=ks.test(msdata\$V4,scrmdata\$V4)
pdf(paste($nsam,\"sampleBLmut\",${t},\"recomb\",${r},\"-KStest.pdf\",sep=\"\"));
par(mfrow=c(2,2))
plot(ecdf(msdata\$V1), xlim=range(c(msdata\$V1, scrmdata\$V1)),col=\"red\", main=\"1st Moment\")
plot(ecdf(scrmdata\$V1), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",m1test\$statistic,sep=\"\"), paste(\"p-value = \",format(m1test\$p.value,scientific = TRUE),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
####
plot(ecdf(msdata\$V2), xlim=range(c(msdata\$V2, scrmdata\$V2)),col=\"red\", main=\"2nd Moment\")
plot(ecdf(scrmdata\$V2), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",m2test\$statistic,sep=\"\"), paste(\"p-value = \",format(m2test\$p.value,scientific = TRUE),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
####
plot(ecdf(msdata\$V3), xlim=range(c(msdata\$V3, scrmdata\$V3)),col=\"red\", main=\"3rd Moment\")
plot(ecdf(scrmdata\$V3), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",m3test\$statistic,sep=\"\"), paste(\"p-value = \",format(m3test\$p.value,scientific = TRUE),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
####
plot(ecdf(msdata\$V4), xlim=range(c(msdata\$V4, scrmdata\$V4)),col=\"red\", main=\"4th Moment\")
plot(ecdf(scrmdata\$V4), add=TRUE, lty=\"dashed\", col=\"blue\")
legend(\"bottomright\",c(paste(\"Tests Statistics = \",m4test\$statistic,sep=\"\"), paste(\"p-value = \",format(m4test\$p.value,scientific = TRUE),sep=\"\")))
legend(\"topleft\",c(\"ms\",\"scrm\"), col=c(\"red\",\"blue\"), pch=16)
####
dev.off();
cat(paste(${nsam},${t},${r},\"|\",
format(m1test\$statistic,digits=4),format(m1test\$p.value,scientific = TRUE),\"|\",
format(m2test\$statistic,digits=4),format(m2test\$p.value,scientific = TRUE),\"|\",
format(m3test\$statistic,digits=4),format(m3test\$p.value,scientific = TRUE),\"|\",
format(m4test\$statistic,digits=4),format(m4test\$p.value,scientific = TRUE),
sep=\"\t\"),file=\"${compareBL}\",append=TRUE);cat(\"\n\",file=\"${compareBL}\",append=TRUE);
" > dummy.r
			R CMD BATCH dummy.r
			done
		done
	done
	
	

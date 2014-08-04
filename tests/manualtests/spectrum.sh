#!/bin/bash
## compare the Spectrum of the segregating sites, frequencies of observing a k mutation...

###########################
#!!!!! need to check is there a formula for the theoretical probability????

mkdir test-spectrum
cd test-spectrum

mst=(10 20 50 100 10)

msNsample=(3 7 10)

msNsample=(3)

rep=1000
#npop=20000

compareSPEC=compareSPEC

rm ${compareSPEC}
#echo -e "compare number of mutations for ${rep} replicates \n\t\t|\t\t\t|\t\t ms \t\t|\t\t scrm\nNsam\ttheta\t|\tmean\tstdv\t|\tmean\tstdv\tstd err\t|\tmean\tstdv \tstd err" >${compareSEG}

for t in "${mst[@]}"
	do
	for nsam in "${msNsample[@]}"
		do
		prefix=${nsam}sample${t}mut
		out=${prefix}out
		nseg=${prefix}Seg
		rm mscolsum*  
		find . -name "xx*" -print0 | xargs -0 rm
		
		ms ${nsam} ${rep} -t ${t} | tail -n +4 | sed '/segsites/d' | sed '/positions/d' | gawk '/^\/\//{f="xx"++d} f{print > f} '
		
		for file in $(seq 1 1 ${rep})
				do 
				sed '/\/\//d' xx${file} | sed 's/.\{1\}/& /g' | awk '
{ for(i=1;i<=NF;++i) t[i]+=$i
  if(n<NF)n=NF
} 
END {
  printf t[1]
  for(i=2;i<=n;++i) printf ","t[i]
  printf "\n"
}
' >> mscolsumsOld
				done
		cat mscolsumsOld | tr '\n' ',' > mscolsums
	
		rm scrmcolsum*  
		find . -name "yy*" -print0 | xargs -0 rm
	
		scrm ${nsam} ${rep} -t ${t} | tail -n +4 | sed '/segsites/d' | sed '/positions/d' | gawk '/^\/\//{f="yy"++d} f{print > f} '
		
        for file in $(seq 1 1 ${rep})
				do 
				sed '/\/\//d' yy${file} | sed 's/.\{1\}/& /g' | awk '
{ for(i=1;i<=NF;++i) t[i]+=$i
  if(n<NF)n=NF
} 
END {
  printf t[1]
  for(i=2;i<=n;++i) printf ","t[i]
  printf "\n"
}
' >> scrmcolsumsOld
				done
		cat scrmcolsumsOld | tr '\n' ',' > scrmcolsums


#echo  -e "Sample size = ${nsam}, theta = ${t}" >> ${compareSPEC}

		
		echo "rm(list=ls());
		#source(\"fun_src.r\");
		msdata=as.numeric(read.table(\"mscolsums\",sep=\",\"));
		msnum=length(msdata)-1
		mstable=table(msdata)
		
		a=as.numeric(read.table(\"scrmcolsums\",sep=\",\"));
		scrmdata=a[!is.na(a)]
		scrmdata=a[a>0]
		scrmnum=length(scrmdata)-1
		scrmtable=table(scrmdata)
		test=chisq.test(cbind(mstable, scrmtable));
		cat(paste(\"Sample size = ${nsam}, theta = ${t}, test statistics = \",
format(test\$statistic,digits=4),\", p-value = \",format(test\$p.value,scientific = TRUE),
sep=\"\"),file=\"${compareSPEC}\",append=TRUE);cat(\"\n\",file=\"${compareSPEC}\",append=TRUE);

		for (i in (1:(${nsam}-1))){
		cat(paste( paste(\"P(S=\",i,\")  \",sep=\"\") , \"|\", format(mstable[i]/msnum,scientific = TRUE), \"|\", format(scrmtable[i]/scrmnum,scientific = TRUE),
sep=\"\t\"),file=\"${compareSPEC}\",append=TRUE);
cat(\"\n\",file=\"${compareSPEC}\",append=TRUE);						
		}
		" > dummy.r
		R CMD BATCH dummy.r
        find . -name "xx*" -print0 | xargs -0 rm
        find . -name "yy*" -print0 | xargs -0 rm
		#rm ms${out} ms${nseg} scrm${out} scrm${nseg}
		done
	done




# make the transpose?
#awk '
#{ 
    #for (i=1; i<=NF; i++)  {
        #a[NR,i] = $i
    #}
#}
#NF>p { p = NF }
#END {    
    #for(j=1; j<=p; j++) {
        #str=a[1,j]
        #for(i=2; i<=NR; i++){
            #str=str" "a[i,j];
        #}
        #print str
    #}
#}' > xx${file}transpose

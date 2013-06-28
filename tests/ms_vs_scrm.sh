#!/bin/bash
msNsample=(2 3 4 7 10 20)

rep=100

## compare TMRCA
compareTMRCA=compareTMRCA

for nsam in "${msNsample[@]}"
	do
	prefix=${nsam}sample
	out=${prefix}out
	tmrca=${prefix}tmrca
	Trees=${prefix}Trees
	ms ${nsam} ${rep} -T | tail -n +4 | grep -v "//" > ms${out}
	cat ms${out} | grep ";" | sed -e 's/\[.*\]//g' > ms${Trees}

	hybrid-Lambda -gt ms${Trees} -tmrca ms${tmrca}


	scrm ${nsam} ${rep} -npop ${npop} -tmrca scrm${tmrca} 
#	echo "rm(list=ls());cat(paste(${nsam},,),file=\"${compareTMRCA}\",append=TRUE);cat(\"\n\",file=\"${compareTMRCA}\",append=TRUE);" > dummy.r
#	R CMD BATCH dummy.r
#	rm ms${out} ms${Trees} ms${tmrca} scrm${tmrca} 
	done


###########################

mst=(10 20 50 100)
msr=(10 20 10 50)
msNsample=(2 3 7 10)

npop=20000
seqlen=100000


	
## compare number of segregating sites




## compare times of genealogy changes



	
	
##cat ${msout} | grep ";" > ${msTrees}
##ntree=$(wc -l ${msTrees} | sed -e "s/${msTrees}//" )
##cut -f1 -d"," ${msTrees} | sed -e 's/\[//g' | sed -e 's/\].*\:/ /g' > ${tmrca}

#echo "rm(list=ls());seqdata=strsplit(scan(\"${segseqfile}\", what=\"charactar\"),\"\");position=scan(\"${position}\")*${seqlen};unlink(\"${VCF}\");" > RgenerateVCF.r
#echo "cat(\"##fileformat=VCFv4.1\n\",file=\"${VCF}\",append=TRUE); " >> RgenerateVCF.r
#echo "cat(\"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\",file=\"${VCF}\",append=TRUE);" >> RgenerateVCF.r
#echo "if (${Nsample}>1){for (i in (1:(${Nsample}-1))){ cat (paste(\"\tNA0000\",i+1,sep=\"\"),file=\"${VCF}\",append=TRUE);}}" >> RgenerateVCF.r
##echo "for (i in (1:(${Nsample}-1))){ cat (paste(\"\tNA0000\",i+1,sep=\"\"),file=\"${VCF}\",append=TRUE);}" >> RgenerateVCF.r
#echo "cat(\"\n\",file=\"${VCF}\",append=TRUE);" >> RgenerateVCF.r
##echo "cat("20	17	.	T	A	3	PASS	NS=3;DP=11;AF=0.017	GT:GQ:DP:HQ	1|0:49:3:58,50\n",file="${VCF}",append=TRUE); " >> RgenerateVCF.r
#echo "for (i in 1:length(position)){cat(1,position[i],\"rs6040355\",\"A\",\"T\",\"67\",\"PASS\",\"NS=2;\",\"GT\",paste(seqdata[[1]][i],\"|\",seqdata[[2]][i],sep=\"\"),sep=\"\t\",file=\"${VCF}\",append=TRUE);" >> RgenerateVCF.r
#echo "if (${Nsample}>1){for (j in (1:(${Nsample}-1))){ cat (paste(\"\t\",seqdata[[2*j+1]][i],\"|\",seqdata[[2*j+2]][i],sep=\"\"),file=\"${VCF}\",append=TRUE);}}
#cat(\"\n\",file=\"${VCF}\",append=TRUE);}" >> RgenerateVCF.r

##echo "for (i in 1:length(position)){cat(1,position[i],\"rs6040355\",\"A\",\"T\",\"67\",\"PASS\",\"NS=2;\",\"GT\",paste(seqdata[[1]][i],\"|\",seqdata[[2]][i],\"\n\",sep=\"\"),sep=\"\t\",file=\"${VCF}\",append=TRUE);}" >> RgenerateVCF.r
##echo "for (i in 1:length(position)){cat(1,position[i],\"rs6040355\",\"A\",\"T\",\"67\",\"PASS\",\"NS=2;\",\"GT:GQ\",paste(seqdata[[1]][i],\"|\",seqdata[[2]][i],\":1,40\n\",sep=\"\"),sep=\"\t\",file=\"${VCF}\",append=TRUE);}" >> RgenerateVCF.r
#echo "pdf(\"${tmrca}.pdf\");plot(c(0,cumsum(read.table(\"${msTreesChangeAt}\")\$V1)), c(read.table(\"${tmrca}\")\$V1*4*${npop},0), xlab=\"site\",ylab=\"${tmrca} / coal unit\",type=\"s\");lines(position,rep(min(read.table(\"${tmrca}\")$V2*4*${npop})*0.9,length(position)),type=\"p\",col=\"red\");dev.off()" >> RgenerateVCF.r

#R CMD BATCH RgenerateVCF.r




#echo "rm(list=ls());TMRCA_file= \"my${TMRCA}\";tmrca=\"${tmrca}\";
#mydata=read.table(TMRCA_file);TMRCA=as.matrix(mydata[,-1]);base=mydata[,1];x=base;numof_y=25;TMRCA_max=max(read.table(tmrca)$V1*4*${npop});
#z=matrix(0, length(base),numof_y);spacing = TMRCA_max/numof_y;y=seq(spacing,TMRCA_max, spacing);for (i in 1:(length(x))){ dummymax=TMRCA_max;
#remaining_num=sum(TMRCA[i,]<dummymax);zindx=numof_y;while(dummymax>0){dummymax=dummymax-spacing; remaining_num2=sum(TMRCA[i,]<dummymax);z[i,zindx]=remaining_num-remaining_num2;zindx=zindx-1;remaining_num=remaining_num2;} };z=z/${Nparticles}; xlim = range(x, finite = TRUE); ylim = range(y, finite = TRUE); zlim = range(z, finite = TRUE); nlevels = 11;levels = pretty(zlim, nlevels);col=topo.colors(length(levels)-1); asp = NA; xaxs = \"i\"; yaxs = \"i\"; las = 1; axes = TRUE; frame.plot = axes; pdf(\"${prefix}combined.pdf\",width=18,height=7); mar.orig <- (par.orig <- par(c(\"mar\", \"las\", \"mfrow\")))\$mar; on.exit(par(par.orig));w <- (3 + mar.orig[2L]) * par(\"csi\") * 2.54;layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w))); par(las = las); mar <- mar.orig; mar[4L] <- mar[2L]; mar[2L] <- 1; par(mar = mar); plot.new(); plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = \"i\", yaxs = \"i\");rect(0, levels[-length(levels)], 1, levels[-1L], col = col);axis(4); mar <- mar.orig; mar[4L] <- 1; par(mar = mar); plot.new();plot.window(xlim, ylim, \"\", xaxs = xaxs, yaxs = yaxs, asp = asp); .filled.contour(x, y, z, levels, col);lines(c(0,cumsum(read.table(\"${msTreesChangeAt}\")\$V1)), c(read.table(tmrca)\$V1*4*${npop},0), xlab=\"site\",ylab=\"tmrca / coal unit\",type=\"s\");title(main = \"\", xlab = \"\", ylab = \"\"); Axis(x, side = 1);Axis(y, side = 2);position=scan(\"${position}\")*${seqlen};lines(position,rep(TMRCA_max*0.9,length(position)),type=\"p\",col=\"red\");dev.off()" > combiningtmrcas.r


#prefix=${nsam}sample${t}mutation${r}recomb
#tmrca=${prefix}tmrca
#Trees=${prefix}Trees
#TreesLast=${prefix}TreesLast
#TMRCA=${prefix}TMRCA
#segseqfile=${prefix}segseq
#position=${prefix}position
#out=${prefix}out

#ms ${nsam} ${rep} -T -t ${t} -r ${r} ${seqlen} | tail -n +4 | grep -v "//" > ms${out}
##cat ${msout} | tail -${msNsample} > ${segseqfile}
#cat ms${out} | grep "segsites" | sed -e "s/segsites: //" > ms${nseg}
##cat ${msout} | grep "positions" | sed -e "s/positions: //" > ${position}
#cat ms${out} | grep ";" | sed -e 's/\[.*\]//g' > ${msTrees}
#cat ${msout} | grep ";" | sed -e 's/\[//g' | sed -e 's/\].*;//g' > ${msTreesLast}
#ntree=$(wc -l ${msTrees} | sed -e "s/${msTrees}//" )


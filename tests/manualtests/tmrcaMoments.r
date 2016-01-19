rm(list=ls())
args=(commandArgs(TRUE))
print(args)
rep = as.numeric(args[1])
seqlen = as.numeric(args[2])
program = args[3]

for (i in (1:rep) ){
  p=read.table(paste("xx",i,"Trees_freq",sep=""), header = F)$V1/seqlen;
  tmrca=read.table(paste("xx",i,"Treestmrca",sep=""), header = F)$V1
  bl=read.table(paste("xx",i,"Treesbl",sep=""), header = F)$V1

  cat(paste(sum(tmrca*p),sum(tmrca^2*p),sum(tmrca^3*p),sum(tmrca^4*p),sep="\t") ,file=paste(program,"tmrca_moments",sep=""),append=TRUE)
  cat("\n",file=paste(program,"tmrca_moments",sep=""),append=TRUE)
  cat(paste(sum(bl*p),sum(bl^2*p),sum(bl^3*p),sum(bl^4*p),sep="\t") ,file=paste(program,"bl_moments",sep=""),append=TRUE)
  cat("\n",file=paste(program,"bl_moments",sep=""),append=TRUE)
}

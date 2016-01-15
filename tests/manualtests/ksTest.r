rm(list=ls())
args=(commandArgs(TRUE))
print(args)

figuretitle = args[1] #scan("figuretitle",what="");
currentcase = args[2] #scan("current_case",what="");

msdata=read.table(args[3])$V1
scrmdata=read.table(args[4])$V1

test=ks.test(msdata,scrmdata)

pdf(paste(currentcase,figuretitle,".pdf",sep=""))
plot(ecdf(msdata), xlim=range(c(msdata, scrmdata)),col="red", main=paste(currentcase, figuretitle))
plot(ecdf(scrmdata), add=TRUE, lty="dashed", col="blue")
legend("bottomright",c(paste("Tests Statistics = ",test$statistic,sep=""), paste("p-value = ",format(test$p.value,digits=4),sep="")))
legend("topleft",c("ms","scrm"), col=c("red","blue"), pch=16)
dev.off();

keyword = "Not "
if ( test$p.value > 0.05 ){
    keyword = ""
}

print( test$p.value )
print ( paste ( "Status ", keyword, "Ok", sep = "") )


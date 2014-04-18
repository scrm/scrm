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

library(binom)

#tab_file = "MIRET.tab"
#in_file = paste(path,"/",tab_file, sep='')
#out_file = paste(path,"/",tab_file, sep='')
get_Freq_Test_Pval=function(tab){
	if(is.na(tab[1]) | (tab[2]==0 & tab[3]==0)){return (NA)}
	else{ return (binom.test(x=tab[2], n=tab[3], p=tab[1], alternative="greater")$p.value)}
}


RF=read.table(in_file, h=T)
if(method == "IES"){
	REF_upper=binom.confint(x=RF$CTL_IES, n=RF$CTL_MAC+RF$CTL_IES, methods = "exact", conf.level = 0.95)$upper
	binomtable=cbind(REF_upper, RF$CUR_IES, RF$CUR_IES+RF$CUR_MAC)
	pvalues=apply(binomtable, 1, 'get_Freq_Test_Pval')
	RF$padj=p.adjust(pvalues, method="BH" )
	RF$SIGNIFICANT = RF$padj < 0.05 & !is.na(RF$padj) & (RF$CTL_MAC+RF$CTL_IES>0) & (RF$CUR_IES+RF$CUR_MAC>0)
}
if(method == "Boundaries"){
	REF_upper=binom.confint(RF$CTL_LEFT, RF$CTL_MAC+RF$CTL_LEFT, methods = "exact", conf.level = 0.95)$upper
	binomtable=cbind(REF_upper, RF$CUR_LEFT, RF$CUR_LEFT+RF$CUR_MAC)
	pvaluesL=apply(binomtable, 1, 'get_Freq_Test_Pval')
	RF$padj_left=p.adjust(pvaluesL, method="BH" )
	RF$SIGNIFICANT_LEFT = RF$padj_left < 0.05 & !is.na(RF$padj_left) & (RF$CTL_MAC+RF$CTL_LEFT>0) & (RF$CUR_LEFT+RF$CUR_MAC>0)

	REF_upper=binom.confint(RF$CTL_RIGHT, RF$CTL_MAC+RF$CTL_RIGHT, methods = "exact", conf.level = 0.95)$upper
	binomtable=cbind(REF_upper, RF$CUR_RIGHT, RF$CUR_RIGHT+RF$CUR_MAC)
	pvaluesR=apply(binomtable, 1, 'get_Freq_Test_Pval')
	RF$padj_right=p.adjust(pvaluesR, method="BH" )
	RF$SIGNIFICANT_RIGHT = RF$padj_right < 0.05 & !is.na(RF$padj_right) & (RF$CTL_MAC+RF$CTL_RIGHT>0) & (RF$CUR_RIGHT+RF$CUR_MAC>0)
}
write.table(x=RF, file=out_file, col.names=T, row.names=F, quote=F, sep="\t")



#!/usr/bin/env Rscript
####################################
#############PARAMETERS#############
####################################
annotation.file=""#two columns tsv file (ensembl ID and gene symbol)
wd=""#with raw count files from HTseq
conditions=c(rep('high',times=),rep('low',times=))
comparison=c("condition","low","high")
prefix='CDX2_high_low'
Number.Samples=
Number.Conditions=
Qt.CutOff=0.85
minReplicates=
pAdjustMethod=""
GLM=TRUE
####################################
####################################
suppressMessages(library('DESeq2'))
suppressMessages(library("gplots"))
suppressMessages(library("ggplot2"))
AddGeneName<-function(l){
	a<-gDE$symbol[gDE$id %in% l[1]];
	if(length(a) == 0){
		return("NA")
	}else{
		(return(as.character(a[1])))
	}
}
####################################
####################################
sampleTable<-as.data.frame(cbind('sampleName'=c(as.vector(sapply(as.vector(grep(".count",list.files(wd),value=TRUE)),function(x) unlist(strsplit(x,"[.]"))[1]))),'fileName'=grep(".count",list.files(wd),value=TRUE),'condition'=conditions))
print(sampleTable)
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=wd,design = ~condition)
if(GLM == TRUE){
	ddsHTSeq<-DESeq(ddsHTSeq)
}else{
	cat("No GLM\n")
	ddsHTSeq<-DESeq(ddsHTSeq,betaPrior=FALSE)
}
cat("Size factors\n")
print(sizeFactors(ddsHTSeq))
normCount<-counts(ddsHTSeq,normalized=TRUE)
colnames(normCount)<-paste(colnames(normCount),'_nCOUNT',sep="")
rawCount<-counts(ddsHTSeq)
colnames(rawCount)<-paste(colnames(rawCount),'_rCOUNT',sep="")
vsd<-assay(varianceStabilizingTransformation(ddsHTSeq))
colnames(vsd)<-paste(colnames(vsd),'_vsd',sep="")
rld<-assay(rlog(ddsHTSeq))
colnames(rld)<-paste(colnames(rld),'_rld',sep="")
diff.exp<-results(ddsHTSeq,contrast=comparison,pAdjustMethod=pAdjustMethod,cooksCutoff=FALSE)
flag<-c()
if(minReplicates >= 3){
	cat("Flaging outliers\n")
	ddsHTSeq.outliers<-replaceOutliers(ddsHTSeq, minReplicates = minReplicates,cooksCutoff=qf(Qt.CutOff,Number.Conditions,Number.Samples-Number.Conditions))
	flag<-as.data.frame(mcols(ddsHTSeq.outliers)$replace)
	rownames(flag)<-rownames(normCount)
	colnames(flag)<-'FLAG'
}else{
	flag<-as.data.frame(rep(NA,times=nrow(diff.exp)))
	rownames(flag)<-rownames(normCount)
	colnames(flag)<-'FLAG'
}
merged.df<-merge(diff.exp,normCount,by="row.names",all=TRUE)
rownames(merged.df)<-merged.df[,1]
merged.df<-merged.df[,-1]
merged.df<-merge(merged.df,vsd,by="row.names",all=TRUE)
rownames(merged.df)<-merged.df[,1]
merged.df<-merged.df[,-1]
merged.df<-merge(merged.df,rld,by="row.names",all=TRUE)
rownames(merged.df)<-merged.df[,1]
merged.df<-merged.df[,-1]
merged.df<-merge(merged.df,rawCount,by="row.names",all=TRUE)
rownames(merged.df)<-merged.df[,1]
merged.df<-merged.df[,-1]
merged.df<-merge(merged.df,flag,by="row.names",all=TRUE)
colnames(merged.df)[1]<-'GeneID'
if(file.exists(annotation.file)){
	cat("Adding Gene Names\n")
	gDE<-read.csv(annotation.file,header=TRUE,sep="\t")
	colnames(gDE)<-c('GeneID','symbol','name')
	merged.df<-merge(merged.df,gDE[,c(1,2)],by="GeneID",all.x=TRUE)
	#merged.df<-cbind(merged.df,apply(merged.df,1,function(l) AddGeneName(l)))
}
merged.df<-merged.df[order(merged.df$padj,decreasing=FALSE),]
cat("Significant DEG: ",nrow(subset(merged.df,merged.df$padj <= 0.05)),"\n")
write.table(merged.df,file=paste(wd,"/",prefix,"_All_MisRegulated_Genes.csv",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
dists<-dist(t(as.matrix(vsd)))
tm<-as.matrix(dists)
colnames(tm)<-sampleTable$sampleName
rownames(tm)<-sampleTable$sampleName
pdf(paste(wd,"/",prefix,"_Expression_Heatmap.pdf", sep = ""))
heatmap.2(tm,margins=c(15,15),trace="none",density.info="none")
dev.off()
save.image(file=paste(wd, "/env.RData", sep = ""),safe=TRUE)
cat("Done\n")

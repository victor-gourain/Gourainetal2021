suppressMessages(library("edgeR"))
human<-read.csv("",header=TRUE,sep="\t")#three columns tsv file (chr\tstart\tend)
human<-cbind(human,'length'=as.vector(human[,3])-as.vector(human[,2]))
human<-human[,c(1,4)]
human$length<-human$length/1000
colnames(human)[1]<-'gID'
human.raw<-human
sample.size<-c()
for(i in list.files(path="",pattern="*.count")){#directory with count file outputted by HTSeq
	print(i)
	ID<-as.character(i)
	tmp<-read.csv(i,header=FALSE,sep="\t")
	colnames(tmp)<-c('gID','rawcount')
	human.raw<-merge(human.raw,tmp,by='gID',all.x=TRUE)
	sample.size<-c(sample.size,sum(as.vector(tmp$rawcount)))
	names(sample.size)[length(sample.size)]<-ID
	colnames(human.raw)[ncol(human.raw)]<-ID
	pmf<-sum(as.vector(tmp$rawcount))/1000000
	tmp$rawcount<-tmp$rawcount/pmf
	colnames(tmp)<-c('gID',ID)
	human<-merge(human,tmp,by='gID',all.x=TRUE)	
}
human.FPKM<-t(rbind(apply(human,1,function(x) as.numeric(x[-c(1,2)])/as.numeric(x[2]))))
colnames(human.FPKM)<-colnames(human[,-c(1,2)])
rownames(human.FPKM)<-as.vector(human$gID)
##RLE normalization to compare samples
COAD.normfactor<-calcNormFactors(na.omit(human.raw[,-c(1,2)]),sample.size,method="RLE")
human.FPKMn<-human.FPKM
for(i in names(COAD.normfactor)){
	human.FPKMn[,i]<-as.vector(human.FPKMn[,i]/COAD.normfactor[i])
}
human.FPKMn<-na.omit(human.FPKMn)
CDX2_FPKM<-human.FPKMn['ENSG00000165556',]
library(reshape)
CDX2.FPKM<-melt(CDX2_FPKM)
CDX2.FPKM<-cbind(CDX2.FPKM,'gene'=rep('CDX2',times=nrow(CDX2.FPKM)))
library(ggplot2)
ggplot(CDX2.FPKM,aes(x=gene,y=value))+geom_boxplot()
save.image('TCGA_preanalysis.RData')

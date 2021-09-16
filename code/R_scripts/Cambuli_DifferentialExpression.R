###
suppressMessages(library('DESeq2'))
#
ES_E14<-read.csv("ES_E14.count",header=FALSE,sep="\t")
ES_E14<-ES_E14[grep('ENSMUS',ES_E14$V1),]
colnames(ES_E14)<-c('ensGeneID','ES_E14')
##
ES_J1<-read.csv("ES_J1.count",header=FALSE,sep="\t")
ES_J1<-ES_J1[grep('ENSMUS',ES_J1$V1),]
colnames(ES_J1)<-c('ensGeneID','ES_J1')
##
iCDX2_ES_MEF<-read.csv("iCdx2_ES_MEF.count",header=FALSE,sep="\t")
iCDX2_ES_MEF<-iCDX2_ES_MEF[grep('ENSMUS',iCDX2_ES_MEF$V1),]
colnames(iCDX2_ES_MEF)<-c('ensGeneID','iCDX2_ES_MEF')
##
iCDX2_ES_TC<-read.csv("iCdx2_ES_TC_plastic.count",header=FALSE,sep="\t")
iCDX2_ES_TC<-iCDX2_ES_TC[grep('ENSMUS',iCDX2_ES_TC$V1),]
colnames(iCDX2_ES_TC)<-c('ensGeneID','iCDX2_ES_TC')
##
TS_EGFP<-read.csv("TS_EGFP.count",header=FALSE,sep="\t")
TS_EGFP<-TS_EGFP[grep('ENSMUS',TS_EGFP$V1),]
colnames(TS_EGFP)<-c('ensGeneID','TS_EGFP')
##
TS_Rs26<-read.csv("TS_Rs26.count",header=FALSE,sep="\t")
TS_Rs26<-TS_Rs26[grep('ENSMUS',TS_Rs26$V1),]
colnames(TS_Rs26)<-c('ensGeneID','TS_Rs26')
##iCDX2-MEF vs ES-E14
ESCdx2MEF_cambuli.data<-merge(ES_E14,iCDX2_ES_MEF,by='ensGeneID')
rownames(ESCdx2MEF_cambuli.data)<-ESCdx2MEF_cambuli.data$ensGeneID
ESCdx2MEF_cambuli.data<-ESCdx2MEF_cambuli.data[,-1]
ESCdx2MEF_info.df<-data.frame('condition'=c('ES_E14','iCDX2_ES_MEF'),'type'=c('SR','SR'))
rownames(ESCdx2MEF_info.df)<-colnames(ESCdx2MEF_cambuli.data)[c(1,2)]
ESCdx2MEF_dds<-DESeqDataSetFromMatrix(countData=ESCdx2MEF_cambuli.data,colData=ESCdx2MEF_info.df,design = ~ condition)
ESCdx2MEF_rld<-rlogTransformation(ESCdx2MEF_dds)
ESCdx2MEF_res<-data.frame(assay(ESCdx2MEF_rld),avgLogExpr=(assay(ESCdx2MEF_rld)[,2]+assay(ESCdx2MEF_rld)[,1])/2,rLogFC=assay(ESCdx2MEF_rld)[,2]-assay(ESCdx2MEF_rld)[,1])
ESCdx2MEF_res<-ESCdx2MEF_res[order(ESCdx2MEF_res$rLogFC,decreasing=TRUE),]
colnames(ESCdx2MEF_res)<-paste('ESCdx2MEF',colnames(ESCdx2MEF_res),sep="_")
#iCDX2-TC vs ES-J1
ESCdx2TC_cambuli.data<-merge(ES_J1,iCDX2_ES_TC,by='ensGeneID')
rownames(ESCdx2TC_cambuli.data)<-ESCdx2TC_cambuli.data$ensGeneID
ESCdx2TC_cambuli.data<-ESCdx2TC_cambuli.data[,-1]
ESCdx2TC_info.df<-data.frame('condition'=c('ES_J1','iCDX2_ES_TC'),'type'=c('SR','SR'))
rownames(ESCdx2TC_info.df)<-colnames(ESCdx2TC_cambuli.data)[c(1,2)]
ESCdx2TC_dds<-DESeqDataSetFromMatrix(countData=ESCdx2TC_cambuli.data,colData=ESCdx2TC_info.df,design = ~ condition)
ESCdx2TC_rld<-rlogTransformation(ESCdx2TC_dds)
ESCdx2TC_res<-data.frame(assay(ESCdx2TC_rld),avgLogExpr=(assay(ESCdx2TC_rld)[,2]+assay(ESCdx2TC_rld)[,1])/2,rLogFC=assay(ESCdx2TC_rld)[,2]-assay(ESCdx2TC_rld)[,1])
ESCdx2TC_res<-ESCdx2TC_res[order(ESCdx2TC_res$rLogFC,decreasing=TRUE),]
colnames(ESCdx2TC_res)<-paste('ESCdx2TC',colnames(ESCdx2TC_res),sep="_")
#TS-eGFP vs ES-J1
TS_cambuli.data<-merge(TS_EGFP,ES_J1,by='ensGeneID')
rownames(TS_cambuli.data)<-TS_cambuli.data$ensGeneID
TS_cambuli.data<-TS_cambuli.data[,-1]
TS_info.df<-data.frame('condition'=c('ES_J1','TS_EGFP'),'type'=c('SR','SR'))
rownames(TS_info.df)<-colnames(TS_cambuli.data)[c(1,2)]
TS_dds<-DESeqDataSetFromMatrix(countData=TS_cambuli.data,colData=TS_info.df,design = ~ condition)
TS_rld<-rlogTransformation(TS_dds)
TS_res<-data.frame(assay(TS_rld),avgLogExpr=(assay(TS_rld)[,2]+assay(TS_rld)[,1])/2,rLogFC=assay(TS_rld)[,2]-assay(TS_rld)[,1])
TS_res<-TS_res[order(TS_res$rLogFC,decreasing=TRUE),]
colnames(TS_res)<-paste('TS',colnames(TS_res),sep="_")
cambuli.res<-merge(ESCdx2MEF_res,ESCdx2TC_res,by="row.names")
rownames(cambuli.res)<-cambuli.res$Row.names
cambuli.res<-cambuli.res[,-1]
cambuli.res<-merge(cambuli.res,TS_res,by="row.names")
source("Annotation_Rfunctions.R")
cambuli.res.anno<-AddAnnotation.FUN(df=cambuli.res,byEnsemblID=TRUE,EnsemblIDc=1,martID="mmusculus_gene_ensembl")
cambuli.res.anno.su<-subset(cambuli.res.anno,
cambuli.res.anno$ESCdx2MEF_ES_E14 >= 0 &
cambuli.res.anno$ESCdx2MEF_iCDX2_ES_MEF >= 0 &
cambuli.res.anno$ESCdx2TC_ES_J1 >= 0 &
cambuli.res.anno$ESCdx2TC_iCDX2_ES_TC >= 0 &
cambuli.res.anno$TS_TS_EGFP >= 0 &
cambuli.res.anno$TS_ES_J1 >= 0
)

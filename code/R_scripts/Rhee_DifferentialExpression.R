###
suppressMessages(library(DESeq2))
#
ES_cells<-read.csv("ES_cells.count",header=FALSE,sep="\t")
ES_cells<-ES_cells[grep('ENSMUS',ES_cells$V1),]
colnames(ES_cells)<-c('ensGeneID','ES')
OECdx2_ES<-read.csv("OECdx2_ES.count",header=FALSE,sep="\t")
OECdx2_ES<-OECdx2_ES[grep('ENSMUS',OECdx2_ES$V1),]
colnames(OECdx2_ES)<-c('ensGeneID','OECdx2')
TS_cells<-read.csv("TS_cells.count",header=FALSE,sep="\t")
TS_cells<-TS_cells[grep('ENSMUS',TS_cells$V1),]
colnames(TS_cells)<-c('ensGeneID','TS')
OECdx2_TS<-read.csv("OECdx2_TS.count",header=FALSE,sep="\t")
OECdx2_TS<-OECdx2_TS[grep('ENSMUS',OECdx2_TS$V1),]
colnames(OECdx2_TS)<-c('ensGeneID','OECdx2')
#ESCdx2 vs ES
ESCdx2_rhee.data<-merge(ES_cells,OECdx2_ES,by='ensGeneID')
rownames(ESCdx2_rhee.data)<-ESCdx2_rhee.data$ensGeneID
ESCdx2_rhee.data<-ESCdx2_rhee.data[,-1]
ESCdx2_info.df<-data.frame('condition'=c('ES','ESCdx2'),'type'=c('SR','SR'))
rownames(ESCdx2_info.df)<-colnames(ESCdx2_rhee.data)[c(1,2)]
ESCdx2_dds<-DESeqDataSetFromMatrix(countData=ESCdx2_rhee.data,colData=ESCdx2_info.df,design = ~ condition)
ESCdx2_rld<-rlogTransformation(ESCdx2_dds)
ESCdx2_res<-data.frame(assay(ESCdx2_rld),avgLogExpr=(assay(ESCdx2_rld)[,2]+assay(ESCdx2_rld)[,1])/2,rLogFC=assay(ESCdx2_rld)[,2]-assay(ESCdx2_rld)[,1])
ESCdx2_res<-ESCdx2_res[order(ESCdx2_res$rLogFC,decreasing=TRUE),]
colnames(ESCdx2_res)<-paste('ESCdx2',colnames(ESCdx2_res),sep="_")
#TSCdx2 vs ES
TSCdx2_rhee.data<-merge(ES_cells,OECdx2_TS,by='ensGeneID')
rownames(TSCdx2_rhee.data)<-TSCdx2_rhee.data$ensGeneID
TSCdx2_rhee.data<-TSCdx2_rhee.data[,-1]
TSCdx2_info.df<-data.frame('condition'=c('ES','TSCdx2'),'type'=c('SR','SR'))
rownames(TSCdx2_info.df)<-colnames(TSCdx2_rhee.data)[c(1,2)]
TSCdx2_dds<-DESeqDataSetFromMatrix(countData=TSCdx2_rhee.data,colData=TSCdx2_info.df,design = ~ condition)
TSCdx2_rld<-rlogTransformation(TSCdx2_dds)
TSCdx2_res<-data.frame(assay(TSCdx2_rld),avgLogExpr=(assay(TSCdx2_rld)[,2]+assay(TSCdx2_rld)[,1])/2,rLogFC=assay(TSCdx2_rld)[,2]-assay(TSCdx2_rld)[,1])
TSCdx2_res<-TSCdx2_res[order(TSCdx2_res$rLogFC,decreasing=TRUE),]
colnames(TSCdx2_res)<-paste('TSCdx2',colnames(TSCdx2_res),sep="_")
#TS vs ES
TS_rhee.data<-merge(ES_cells,TS_cells,by='ensGeneID')
rownames(TS_rhee.data)<-TS_rhee.data$ensGeneID
TS_rhee.data<-TS_rhee.data[,-1]
TS_info.df<-data.frame('condition'=c('ES','TS'),'type'=c('SR','SR'))
rownames(TS_info.df)<-colnames(TS_rhee.data)[c(1,2)]
TS_dds<-DESeqDataSetFromMatrix(countData=TS_rhee.data,colData=TS_info.df,design = ~ condition)
TS_rld<-rlogTransformation(TS_dds)
TS_res<-data.frame(assay(TS_rld),avgLogExpr=(assay(TS_rld)[,2]+assay(TS_rld)[,1])/2,rLogFC=assay(TS_rld)[,2]-assay(TS_rld)[,1])
TS_res<-TS_res[order(TS_res$rLogFC,decreasing=TRUE),]
colnames(TS_res)<-paste('TS',colnames(TS_res),sep="_")
Rhee.res<-merge(ESCdx2_res,TSCdx2_res,by="row.names")
rownames(Rhee.res)<-Rhee.res$Row.names
Rhee.res<-Rhee.res[,-1]
Rhee.res<-merge(Rhee.res,TS_res,by="row.names")
Rhee.res.anno<-AddAnnotation.FUN(df=Rhee.res,byEnsemblID=TRUE,EnsemblIDc=1,martID="mmusculus_gene_ensembl")
Rhee.res.anno.su<-subset(Rhee.res.anno,
Rhee.res.anno$ESCdx2_ES >=0 &
Rhee.res.anno$ESCdx2_OECdx2 >=0 &
Rhee.res.anno$TSCdx2_ES >=0 &
Rhee.res.anno$TSCdx2_OECdx2 >=0 &
Rhee.res.anno$TS_ES >=0 &
Rhee.res.anno$TS_TS >=0
)
############################
############################
#ESCdx2 vs TS
ESCdx2_rhee.data<-merge(TS_cells,OECdx2_ES,by='ensGeneID')
rownames(ESCdx2_rhee.data)<-ESCdx2_rhee.data$ensGeneID
ESCdx2_rhee.data<-ESCdx2_rhee.data[,-1]
ESCdx2_info.df<-data.frame('condition'=c('TS','ESCdx2'),'type'=c('SR','SR'))
rownames(ESCdx2_info.df)<-colnames(ESCdx2_rhee.data)[c(1,2)]
ESCdx2_dds<-DESeqDataSetFromMatrix(countData=ESCdx2_rhee.data,colData=ESCdx2_info.df,design = ~ condition)
ESCdx2_rld<-rlogTransformation(ESCdx2_dds)
ESCdx2_res<-data.frame(assay(ESCdx2_rld),avgLogExpr=(assay(ESCdx2_rld)[,2]+assay(ESCdx2_rld)[,1])/2,rLogFC=assay(ESCdx2_rld)[,2]-assay(ESCdx2_rld)[,1])
ESCdx2_res<-ESCdx2_res[order(ESCdx2_res$rLogFC,decreasing=TRUE),]
colnames(ESCdx2_res)<-paste('ESCdx2',colnames(ESCdx2_res),sep="_")
#TSCdx2 vs TS
TSCdx2_rhee.data<-merge(TS_cells,OECdx2_TS,by='ensGeneID')
rownames(TSCdx2_rhee.data)<-TSCdx2_rhee.data$ensGeneID
TSCdx2_rhee.data<-TSCdx2_rhee.data[,-1]
TSCdx2_info.df<-data.frame('condition'=c('TS','TSCdx2'),'type'=c('SR','SR'))
rownames(TSCdx2_info.df)<-colnames(TSCdx2_rhee.data)[c(1,2)]
TSCdx2_dds<-DESeqDataSetFromMatrix(countData=TSCdx2_rhee.data,colData=TSCdx2_info.df,design = ~ condition)
TSCdx2_rld<-rlogTransformation(TSCdx2_dds)
TSCdx2_res<-data.frame(assay(TSCdx2_rld),avgLogExpr=(assay(TSCdx2_rld)[,2]+assay(TSCdx2_rld)[,1])/2,rLogFC=assay(TSCdx2_rld)[,2]-assay(TSCdx2_rld)[,1])
TSCdx2_res<-TSCdx2_res[order(TSCdx2_res$rLogFC,decreasing=TRUE),]
colnames(TSCdx2_res)<-paste('TSCdx2',colnames(TSCdx2_res),sep="_")

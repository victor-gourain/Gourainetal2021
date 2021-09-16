##R functions
####
##Find overlapp between genomic coordinates
####
suppressMessages(library('GenomicRanges'))
suppressMessages(library('IRanges'))
suppressMessages(library('GenomicFeatures'))
oCT=10
ByOverlap.FUN<-function(
gr.object=NA,
gr.annotation=NA,
st=NA
){
	tmp.0<-as.data.frame(gr.annotation,row.names=NULL)
	gr.df<-as.data.frame(gr.object)
	gr.df<-cbind(gr.df,'id'=paste("tagid_",1:nrow(gr.df),sep=""))
	annotation.df<-as.data.frame(gr.annotation)
	annotation.df<-as.data.frame(findOverlaps(gr.object,gr.annotation,ignore.strand=st,type="any",select="arbitrary",minoverlap=oCT))
	annotation.df<-cbind(annotation.df,'id'=paste("tagid_",1:nrow(annotation.df),sep=""))
	tmp.1<-subset(tmp.0,rownames(tmp.0) %in% annotation.df[,1])
	tmp.2<-cbind(rownames(tmp.1),tmp.1)
	colnames(tmp.2)[1]<-colnames(annotation.df)[1]
	M.0<-merge(annotation.df,tmp.2,by=colnames(annotation.df)[1],all.x=TRUE)
	colnames(M.0)<-c('annoRoW',colnames(M.0)[2],paste('Annotation',colnames(M.0)[-c(1,2)],sep="_"))
	M.1<-merge(gr.df,M.0[,-1],by='id',all.x=TRUE)
	return(M.1[,-1])
}
##Annotate feature directly before genomic coordinates
precede.FUN<-function(
	gr.object=NA,
	gr.annotation=NA
){
	tmp.0<-as.data.frame(gr.annotation,row.names=NULL)	
	gr.df<-as.data.frame(gr.object)
	gr.df<-cbind(gr.df,'id'=paste("tagid_",1:nrow(gr.df),sep=""))
	annotation.df<-as.data.frame(gr.annotation)	
	annotation.df<-as.data.frame(precede(gr.object,gr.annotation,ignore.strand=TRUE))
	annotation.df<-cbind(annotation.df,'id'=paste("tagid_",1:nrow(annotation.df),sep=""))
	tmp.1<-subset(tmp.0,rownames(tmp.0) %in% annotation.df[,1])
	tmp.2<-cbind(rownames(tmp.1),tmp.1)
	colnames(tmp.2)[1]<-colnames(annotation.df)[1]
	M.0<-merge(annotation.df,tmp.2,by=colnames(annotation.df)[1],all.x=TRUE)
	colnames(M.0)<-c('annoRoW',colnames(M.0)[2],paste('Annotation_following',colnames(M.0)[-c(1,2)],sep="_"))
	M.1<-merge(gr.df,M.0[,-1],by='id',all.x=TRUE)
	return(M.1)
}
##follow.FUN: to find the annotation feature righ after an intergenic query, MUST BE STRANDED!!!
follow.FUN<-function(
	gr.object=NA,
	gr.annotation=NA
){
	tmp.0<-as.data.frame(gr.annotation,row.names=NULL)	
	gr.df<-as.data.frame(gr.object)
	gr.df<-cbind(gr.df,'id'=paste("tagid_",1:nrow(gr.df),sep=""))
	annotation.df<-as.data.frame(gr.annotation)	
	annotation.df<-as.data.frame(follow(gr.object,gr.annotation,ignore.strand=TRUE))
	annotation.df<-cbind(annotation.df,'id'=paste("tagid_",1:nrow(annotation.df),sep=""))
	tmp.1<-subset(tmp.0,rownames(tmp.0) %in% annotation.df[,1])
	tmp.2<-cbind(rownames(tmp.1),tmp.1)
	colnames(tmp.2)[1]<-colnames(annotation.df)[1]
	M.0<-merge(annotation.df,tmp.2,by=colnames(annotation.df)[1],all.x=TRUE)
	colnames(M.0)<-c('annoRoW',colnames(M.0)[2],paste('Annotation_preceding',colnames(M.0)[-c(1,2)],sep="_"))
	M.1<-merge(gr.df,M.0[,-1],by='id',all.x=TRUE)
	return(M.1)
}

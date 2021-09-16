version<-'v3'
date<-'07Apr17'
author<-'V.GOURAIN'
routines<-"AddAnnotation.FUN,TSSAnnotation.FUN,dist.FUN,GOEnrichment.FUN,ExternalEnrichment.FUN,ortho.FUN"
cat("Version:",version," Date: ",date," Author: ",author," Functions: ",routines,"\n")
#######################
##Annotation function##
#######################
##Input: a data frame, a boolean about the type of annotation, minimal information, a bioMart ID
##Output: a data frame
##Comments:
AddAnnotation.FUN<-function(
	df=NA,
	byCoordinates=FALSE,
	byEnsemblID=FALSE,	
	EnsemblIDc=NA,
	seqnamec=NA,
	startc=NA,
	endc=NA,
	strandc=NA,
	martID=NA,
	hostID='www.ensembl.org'
){
	#"drerio_gene_ensembl"
	suppressMessages(library('GenomicRanges'))
	suppressMessages(library('IRanges'))
	suppressMessages(library('GenomicFeatures'))
	suppressMessages(library('biomaRt'))		
	ignore.strand=TRUE
	##Test on parameters
	cat("\tAddAnnotation: Tests on parameters.\n")
	if(is.na(martID)){
		stop("\tAddAnnotation ERROR: A mart dataset is required.\n")
	}
	if((! is.matrix(df))&&(! is.data.frame(df))){
		stop("\tAddAnnotation ERROR: The data provided are not data frame or matrix. Exit function.\n")
	}else{
		if(byEnsemblID){
			if(is.na(EnsemblIDc) || ! EnsemblIDc <= ncol(df) || ! all(grepl("ENS[A-Z]+[0-9]{11}",df[,EnsemblIDc],perl=TRUE))){
				stop("\tAddAnnotation ERROR: These are not Ensembl ID like ENS[A-Z]+[0-9]{11} or Ensembl ID column number out of range. Exit function.\n")
			}
		}else{
			if(byCoordinates){
				if(is.na(seqnamec) || is.na(startc)){
					stop("\tAddAnnotation ERROR: The minimum requirement for annotation by coordinates is a seqname and a start position\n")
				}else{
					if(! seqnamec <= ncol(df) || ! all(grepl("chr[0-9A-Za-z]{1,10}",df[,seqnamec],perl=TRUE)) || ! startc <= ncol(df) || ! is.numeric(df[,startc])){
						stop("\tAddAnnotation ERROR: Start position not numeric or seqname out of range or seqnames are not formated.\n")
					}else{
						colnames(df)[seqnamec]<-'seqname'
						colnames(df)[startc]<-'start'
						if(is.na(endc)){
							cat("\tAddAnnotation WARNING: Use startc as endc\n")
							df<-cbind(df,'end'=df[,startc])
							endc<-ncol(df)
						}else{
							if(! seqnamec <= ncol(df) || ! is.numeric(df[,endc])){
								stop("\tAddAnnotation ERROR: End position is not numeric or out of range.\n")
							}else{
								colnames(df)[endc]<-'end'
							}
						}
						if(is.na(strandc)){
							df<-cbind(df,'strand'=rep('*',times=nrow(df)))
							strandc=ncol(df)
							ignore.strand=TRUE
							cat("\tAddAnnotation WARNING: Unstranded annotation.\n")
						}else{	
							if(! strandc <= ncol(df)){
								stop("\tAddAnnotation ERROR: strandc out of range.\n")
							}else{
								if(all(grepl("[+-]",df[,strandc],perl=TRUE))){
									colnames(df)[strandc]<-'strand'
									ignore.strand=FALSE
								}else{
									stop("\tAddAnnotation ERROR: strand has to be either + ot -.\n")
								}
							}
						}
					}
				}
			}else{
				stop("\tAddAnnotation ERROR: Must by eiter byEnsemblID or byCoordinates.\n")
			}	
		}
	}
	##Reformating the annotation data
	cat("\tAddAnnotation: Downloading and formating annotation data.\n")
	ensembl<-useMart('ENSEMBL_MART_ENSEMBL',host=hostID,dataset=martID)
	#OR#ensembl=useMart("ensembl",dataset=martID)	#to use if the previous line rise error
	geneDataExternaltemp<-getBM(attributes=c('ensembl_gene_id','entrezgene_id','external_gene_name','description','chromosome_name','start_position','end_position','strand'), mart=ensembl)
	geneDataExternaltemp[,'chromosome_name']<-paste('chr',geneDataExternaltemp[,'chromosome_name'],sep="")
	geneDataExternaltemp[geneDataExternaltemp$strand == -1,'strand']<-'-'
	geneDataExternaltemp[geneDataExternaltemp$strand == 1,'strand']<-'+'
	geneDataExternaltemp<-cbind('BioMartEntry'=paste('BME',rownames(geneDataExternaltemp),sep="_"),geneDataExternaltemp)	
	colnames(geneDataExternaltemp)<-c('BioMartEntry','ensembl_Gene_ID','RefSeq_Gene_ID','Gene_Symbole','Gene_description','seqname','start','end','strand')
	geneDataExternaltemp[geneDataExternaltemp == ""]<-NA
	geneDataExternaltemp<-subset(geneDataExternaltemp, !is.na(geneDataExternaltemp$ensembl_Gene_ID) || !is.na(geneDataExternaltemp$seqname) || !is.na(geneDataExternaltemp$start) || !is.na(geneDataExternaltemp$end))
	annotation<-geneDataExternaltemp
	##Functions
	Func1<-function(l,EnsemblIDc){
		a<-annotation$Gene_Symbole[annotation$ensembl_Gene_ID %in% l[EnsemblIDc]]
		if(length(a) == 0){
			return("NA")
		}else{
			if(is.vector(a)){
				return(paste(a,collapse="/"))
			}else{
				return(as.character(a))
			}
		}
	}
	transFESSE.FUN<-function(gr.object,gr.annotation,st){
		tmp.0<-as.data.frame(gr.annotation,row.names=NULL)	
		gr.df<-as.data.frame(gr.object)
		gr.df<-cbind(gr.df,'id'=paste("tagid_",1:nrow(gr.df),sep=""))
		annotation.df<-as.data.frame(gr.annotation)	
		annotation.df<-as.data.frame(nearest(gr.object,gr.annotation,ignore.strand=st))
		annotation.df<-cbind(annotation.df,'id'=paste("tagid_",1:nrow(annotation.df),sep=""))
		tmp.1<-subset(tmp.0,rownames(tmp.0) %in% annotation.df[,1])
		tmp.2<-cbind(rownames(tmp.1),tmp.1)
		colnames(tmp.2)[1]<-colnames(annotation.df)[1]
		M.0<-merge(annotation.df,tmp.2,by=colnames(annotation.df)[1],all.x=TRUE)
		colnames(M.0)<-c('annoRoW',colnames(M.0)[2],paste('Annotation',colnames(M.0)[-c(1,2)],sep="_"))
		M.1<-merge(gr.df,M.0[,-1],by='id',all.x=TRUE)
		return(M.1[,-1])
	}
	##Main
	cat("\tAddAnnotation: Annotation.\n")
	if(byCoordinates){
		df.gr<-GRanges(seqnames=df$seqname,ranges=IRanges(start=df$start,end=df$end),strand=df$strand)
		df.values<-as.data.frame(df[,-c(seqnamec,startc,endc,strandc),drop=FALSE])		
		values(df.gr)<-cbind(values(df.gr),df.values)		
		annotation.gr<-GRanges(seqnames=annotation$seqname,ranges=IRanges(start=annotation$start,end=annotation$end),strand=annotation$strand)
		annotation.values<-as.data.frame(annotation[,-c(ncol(annotation)-3,ncol(annotation)-2,ncol(annotation)-1,ncol(annotation)),drop=FALSE])
		values(annotation.gr)<-cbind(values(annotation.gr),annotation.values)
		df.annotated<-transFESSE.FUN(df.gr,annotation.gr,st=ignore.strand)
		return(df.annotated)
	}else{
		df.annotated<-cbind(df,apply(df,1,function(x) Func1(x,EnsemblIDc)))
		colnames(df.annotated)[ncol(df.annotated)]<-'Gene_Symbole'
		return(df.annotated)
	}
}
###########################
##TSS annotation function##
###########################
##Input: data frame of TSS, martID, CT on distance
##Output: a list of two data frames named 'promoters' and 'others'
##Comments:depend on AddAnnotation.FUN, problem with chromosomes
TSSAnnotation.FUN<-function(
	martID=NA,
	df=NA,
	seqnamec=NA,
	startc=NA,
	endc=NA,
	strandc=NA,
	CT=500,
	chromosome.number=NA
){
	suppressMessages(library('GenomicRanges'))
	suppressMessages(library('IRanges'))
	suppressMessages(library('GenomicFeatures'))
	suppressMessages(library('biomaRt'))
	##Test on parameters
	cat("\tTSSAnnotation: Tests on parameters.\n")
	if(is.na(martID)){
		stop("\tTSSAnnotation ERROR: A mart dataset is required.\n")
	}
	if((! is.matrix(df))&&(! is.data.frame(df))){
		stop("\tTSSAnnotation ERROR: The data provided are not data frame or matrix. Exit function.\n")
	}else{
		if(is.na(seqnamec) || is.na(startc)){
			stop("\tTSSAnnotation ERROR: The minimum requirement for annotation by coordinates is a seqname and a start position\n")
		}else{
			if(! seqnamec <= ncol(df) || ! all(grepl("chr[0-9]{1,2}",df[,seqnamec],perl=TRUE)) || ! startc <= ncol(df) || ! is.numeric(df[,startc])){
				stop("\tTSSAnnotation ERROR: Start position not numeric or seqname out of range or seqnames are not formated.\n")
			}else{
				colnames(df)[seqnamec]<-'seqname'
				colnames(df)[startc]<-'start'
				if(is.na(endc)){
					cat("\tTSSAnnotation WARNING: Use startc as endc\n")
					df<-cbind(df,'end'=df[,startc])
					endc<-ncol(df)
				}else{
					if(! seqnamec <= ncol(df) || ! is.numeric(df[,endc])){
						stop("\tTSSAnnotation ERROR: End position is not numeric or out of range.\n")
					}else{
						colnames(df)[endc]<-'end'
					}
				}
				if(is.na(strandc)){
					stop("\tTSSAnnotation ERROR: Unstranded annotation is no supported.\n")
				}else{	
					if(! strandc <= ncol(df)){
						stop("\tTSSAnnotation ERROR: strandc out of range.\n")
					}else{
						if(all(grepl("[+-]",df[,strandc],perl=TRUE))){
							colnames(df)[strandc]<-'strand'
						}else{
							stop("\tTSSAnnotation ERROR: strand has to be either + ot -.\n")
						}
					}
				}
			}
		}
	}
	##Reformating the annotation data
	cat("\tTSSAnnotation: Downloading and formating promoter regions.\n")
	ensembl<-useMart('ENSEMBL_MART_ENSEMBL',host="www.ensembl.org",dataset=martID)
	#ensembl=useMart("ensembl",dataset=martID)
	transcriptDataExternaltemp<-getBM(attributes=c('ensembl_transcript_id','ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position','strand','transcription_start_site'), mart=ensembl)
	transcriptDataExternaltemp[,'chromosome_name']<-paste('chr',transcriptDataExternaltemp[,'chromosome_name'],sep="")
	CHR<-paste('chr',seq(1,chromosome.number),sep="")
	transcriptDataExternaltemp<-subset(transcriptDataExternaltemp,transcriptDataExternaltemp$chromosome_name %in% CHR)
	transcriptDataExternaltemp[transcriptDataExternaltemp$strand == -1,'strand']<-'-'
	transcriptDataExternaltemp[transcriptDataExternaltemp$strand == 1,'strand']<-'+'
	TSSReformating.FUN<-function(
		transcriptID=NA,
		chr=NA,
		strand,
		start=NA,
		end,
		tss
	){
		#if a gene exists with at least a chromosome and a start
		if((!is.na(transcriptID))&&(!is.na(chr))&&(!is.na(start))){
			#if a tss exists, well a tss exists
			if(!is.na(tss)){
				return(tss)
			}
			#if no strand or end but a start, tss is the start
			if(is.na(strand)||is.na(end)){
				return(start)
			}else{
				#if minus strand, tss is end
				if(strand == '-'){
					return(end)
				}else{
					#if plus strand, tss is start
					if(strand == '+'){
						return(start)
					}else{
						return(NA)
					}
				}
			}	
		}else{
			return(NA)
		}
	}
	tssrf<-cbind(transcriptDataExternaltemp,'Putative_TSS'=apply(transcriptDataExternaltemp,1,function(x) TSSReformating.FUN(transcriptID=x[1],chr=x[4],start=as.numeric(x[5]),end=as.numeric(x[6]),strand=x[7],tss=as.numeric(x[8]))))
	promoterRegion.FUN<-function(
		transcriptID=NA,
		geneID=NA,
		geneName=NA,
		chr=NA,
		strand=NA,
		tss=NA
	){
		#if a gene exists with at least ID chr tss and strand
		if((!is.na(geneID))&&(!is.na(chr))&&(!is.na(tss))&&(!is.na(strand))){
			#if tss is numeric
			if((is.numeric(tss))&&(strand == '-' || strand == '+')){
				start<- tss - CT
				end<- tss + CT
				if(start <= 0){
					start<- 1
				}
				if(start < end){
					return(c(transcriptID,geneID,geneName,chr,start,end,strand))
				}
			}
		}
	}
	tss.df<-as.data.frame(t(apply(tssrf,1,function(x) promoterRegion.FUN(transcriptID=x[1],geneID=x[2],geneName=x[3],chr=x[4],tss=as.numeric(x[9]),strand=x[7]))))
	colnames(tss.df)<-c('ensTranscriptID','ensembl_gene_id','external_gene_name','seqnames','start','end','strand')
	tss.df$start<-as.numeric(as.vector(tss.df$start))
	tss.df$end<-as.numeric(as.vector(tss.df$end))
	tss.gr<-GRanges(seqnames=tss.df$seqnames,ranges=IRanges(start=tss.df$start,end=tss.df$end),strand=tss.df$strand)
	tss.values<-as.data.frame(tss.df[,c('ensTranscriptID','ensembl_gene_id','external_gene_name'),drop=FALSE])		
	values(tss.gr)<-cbind(values(tss.gr),tss.values)
	##Reformating the input data
	cat("\tTSSAnnotation: Formating the input data.\n")
	df.gr<-GRanges(seqnames=df$seqname,ranges=IRanges(start=df$start,end=df$end),strand=df$strand)
	df.values<-as.data.frame(df[,-c(seqnamec,startc,endc,strandc),drop=FALSE])		
	values(df.gr)<-cbind(values(df.gr),df.values)
	promoterByOverlap.FUN<-function(
		gr.object=NA,
		gr.annotation=NA
	){
		tmp.0<-as.data.frame(gr.annotation,row.names=NULL)	
		gr.df<-as.data.frame(gr.object)
		gr.df<-cbind(gr.df,'id'=paste("tagid_",1:nrow(gr.df),sep=""))
		annotation.df<-as.data.frame(gr.annotation)	
		annotation.df<-as.data.frame(findOverlaps(gr.object,gr.annotation,type="any",select="arbitrary"))
		annotation.df<-cbind(annotation.df,'id'=paste("tagid_",1:nrow(annotation.df),sep=""))
		tmp.1<-subset(tmp.0,rownames(tmp.0) %in% annotation.df[,1])
		tmp.2<-cbind(rownames(tmp.1),tmp.1)
		colnames(tmp.2)[1]<-colnames(annotation.df)[1]
		M.0<-merge(annotation.df,tmp.2,by=colnames(annotation.df)[1],all.x=TRUE)
		colnames(M.0)<-c('annoRoW',colnames(M.0)[2],paste('Annotation',colnames(M.0)[-c(1,2)],sep="_"))
		M.1<-merge(gr.df,M.0[,-1],by='id',all.x=TRUE)
		return(M.1[,-1])
	}
	##Annotataion
	cat("\tTSSAnnotation: Annotation of TSS.\n")
	tss.annotated<-promoterByOverlap.FUN(gr.object=df.gr,gr.annotation=tss.gr)
	tss.annotated.na<-tss.annotated[!complete.cases(tss.annotated),!grepl('Annotation',colnames(tss.annotated))]
	tss.annotated.promoters<-tss.annotated[complete.cases(tss.annotated),]
	tss.annotated.na<-subset(tss.annotated.na,tss.annotated.na$seqnames %in% CHR)
	tss.annotated.na$seqnames<-as.vector(tss.annotated.na$seqnames)
	tss.non.promoter<-AddAnnotation.FUN(df=tss.annotated.na,byCoordinates=TRUE,seqnamec=1,startc=2,endc=3,strandc=5,martID='drerio_gene_ensembl')
	annotated.list<-list()
	annotated.list[['promoters']]<-tss.annotated.promoters
	annotated.list[['others']]<-tss.non.promoter
	return(annotated.list)
}
#####################
##Distance function##
#####################
##Input: qstart qend astar aend qstrand
##Output: a single value, NA or numeric
##Comments:
dist.FUN<-function(Qstart,Qend,Astart,Aend,Qstrand){
		if(is.na(Qstart) || is.na(Qend) || is.na(Astart) || is.na(Aend)){
			return(NA)
		}else{
			##if either Qstart or Qend is between Astart and Aend it overlaps
			if((Qstart >= Astart && Qstart <= Aend)||(Qend >= Astart && Qend <= Aend)){
				return(0)
			}else{
				##else if strand is +
				if(Qstrand == '+'){
					if(Qend < Astart){
						return(Astart-Qend)
					}else{
						if(Aend < Qstart){
							return(Aend - Qstart)
						}else{
							return(NA)
						}
					}
				}else{
					##else if strand is -
					if(Qend < Astart){
						return(Qend - Astart)
					}else{
						if(Aend < Qstart){
							return(Qstart - Aend)
						}else{
							return(NA)
						}
					}
				}
			}
		}
	}
######################
##GO Term enrichment##
######################
##Input: a vector of ensembl gene ID, a bioMart ID, a gene universe and a multiple testing method 
##Output: a data frame
##Comments:
GOEnrichment.FUN<-function(
	gene.list=NA,
	martID=NA,
	gene.universe.default=NA,
	multiple.testing.method=NA
){
	cat("\tGO Terms Enrichment: Tests on parameters.\n")
	if(is.na(martID)){
		stop("\tGO Terms Enrichment ERROR: A mart dataset is required.\n")
	}
	if(is.na(gene.universe.default)||!is.numeric(gene.universe.default)){
		stop("\tGO Terms Enrichment ERROR: A gene universe is required and should be numeric.\n")
	}
	if(is.na(multiple.testing.method)||!is.element(multiple.testing.method,c('bonferroni','holm','hochberg','hommel','BH','fdr','BY'))){
		stop("\tGO Terms Enrichment ERROR: A multiple testing method is required.\n")
	}
	if(is.na(gene.list)||!is.vector(gene.list)){
		stop("\tGO Terms Enrichment ERROR: A gene list as vector is requiered.\n")	
	}else{
		if(! all(grepl("ENS[A-Z]+[0-9]{11}",gene.list,perl=TRUE))){
	 		stop("\tGO Terms Enrichment ERROR: These are not ensembl gene ID.\n")
		}
	}
	suppressMessages(library(biomaRt))
	cat("\tGO Terms Enrichment: Downloading and formating Gene Ontology data.\n")
	ensembl<-useMart('ENSEMBL_MART_ENSEMBL',host="www.ensembl.org",dataset=martID)
	#ensembl=useMart("ensembl",dataset=martID)	#to use if the previous line rise error
	GODataExternaltemp<-getBM(attributes=c('ensembl_gene_id','external_gene_name','go_id','name_1006'), mart=ensembl)
	GODataExternaltemp<-cbind('BioMartEntry'=paste('BM_GO',rownames(GODataExternaltemp),sep="_"),GODataExternaltemp)	
	colnames(GODataExternaltemp)<-c('GO_BioMartEntry','GO_ensembl_Gene_ID','GO_Gene_Symbole','GO_GO_id','GO_GO_name')
	GODataExternaltemp[GODataExternaltemp == ""]<-NA
	GODataExternaltemp<-subset(GODataExternaltemp, !is.na(GODataExternaltemp$GO_ensembl_Gene_ID) || !is.na(GODataExternaltemp$GO_GO_id))
	GO.data<-GODataExternaltemp
	cat("\tGO Terms Enrichment: Gene universe: ")
	gene.universe<-gene.universe.default
	if(length(unique(GO.data$GO_ensembl_Gene_ID))>=gene.universe.default){
		gene.universe<-length(unique(GO.data$GO_ensembl_Gene_ID))	
		cat(gene.universe,"\n")
	}else{
		gene.universe<-gene.universe.default
		cat(gene.universe,"\n")
	}
	##Number of associated genes to each X in the gene universe
	total.genes.per.GO<-table(GO.data$GO_GO_id)
	total.genes.per.GO<-total.genes.per.GO[total.genes.per.GO != 0]
	##Number of associated genes to each X in the gene list
	GO.data.sub<-subset(GO.data,GO.data$GO_ensembl_Gene_ID %in% gene.list)
	list.genes.per.GO<-table(GO.data.sub$GO_GO_id)
	list.genes.per.GO<-list.genes.per.GO[list.genes.per.GO != 0]
	bigMat<-data.frame(ID=character(),NAME=character(),Fa=numeric(), Fb=numeric(), Fc=numeric(), Fd=numeric(), FPV=numeric(), GENE=character(),stringsAsFactors=FALSE)
	cat("\tGO Terms Enrichment: Statistical analysis.\n")
	for(GO in unique(GO.data.sub$GO_GO_id)){
		tempmat<-matrix(c(list.genes.per.GO[GO],(length(gene.list)-list.genes.per.GO[GO]),total.genes.per.GO[GO],(gene.universe-total.genes.per.GO[GO])),ncol = 2, nrow = 2)
		if((is.na(list.genes.per.GO[GO])==TRUE)||(is.na(length(gene.list)-list.genes.per.GO[GO]) == TRUE)||(is.na(total.genes.per.GO[GO]) == TRUE)||(is.na(gene.universe-total.genes.per.GO[GO]) == TRUE)){
			print(tempmat)
			print(GO)
			print("NA value")
			next
		}
		if((is.numeric(list.genes.per.GO[GO])==FALSE)||(is.numeric(length(gene.list)-list.genes.per.GO[GO]) == FALSE)||(is.numeric(total.genes.per.GO[GO]) == FALSE)||(is.numeric(gene.universe-total.genes.per.GO[GO]) == FALSE)){
			print(tempmat)
			print(GO)
			print("Non numeric value")
			next
		}
		if((is.infinite(list.genes.per.GO[GO])==TRUE)||(is.infinite(length(gene.list)-list.genes.per.GO[GO]) == TRUE)||(is.infinite(total.genes.per.GO[GO]) == TRUE)||(is.infinite(gene.universe-total.genes.per.GO[GO]) == TRUE)){
			print(tempmat)
			print(GO)
			print("Infinite value")
			next
		}
		if((list.genes.per.GO[GO] < 0)||((length(gene.list)-list.genes.per.GO[GO]) < 0)||(total.genes.per.GO[GO] < 0)||((gene.universe-total.genes.per.GO[GO]) < 0)){
			print(tempmat)
			print(GO)
			print("Negative value")
			next
		}
		ft<-fisher.test(tempmat,alternative="greater")
		tempmat[1,1]<-tempmat[1,1]-1
		GO.name<-unique(subset(GO.data,GO.data$GO_GO_id == GO)[,c('GO_GO_id','GO_GO_name')])
		associated.genes<-paste(subset(GO.data.sub,GO.data.sub$GO_GO_id == GO)[,'GO_ensembl_Gene_ID'],subset(GO.data.sub,GO.data.sub$GO_GO_id == GO)[,'GO_Gene_Symbole'],sep=":",collapse=",")
		bigMat<-rbind(bigMat,data.frame(ID=GO,NAME=GO.name,Fa=list.genes.per.GO[GO],Fb=(length(gene.list)-list.genes.per.GO[GO]),Fc=total.genes.per.GO[GO],Fd=(gene.universe-total.genes.per.GO[GO]),FPV=ft$p.value,GENE=associated.genes))
	}
	bigMat.adjusted<-cbind(bigMat,'adjusted_pvalue'=p.adjust(as.vector(bigMat[,'FPV']),method=multiple.testing.method))
	return(bigMat.adjusted)
}
#######################
##External enrichment##
#######################
##Input: a vector of ID, a universe value, a multiple testing method and a file of associations 
##Output: a data frame
##input file format: with header, tab-delimited, ExtID\tExtName\tIntID\tIntName
##Comments:
ExternalEnrichment.FUN<-function(
	gene.list=NA,
	gene.universe.default=NA,
	multiple.testing.method=NA,
	associations.file=NA
){
	cat("\tExternal Enrichment: Tests on parameters.\n")
	if(is.na(gene.universe.default)||!is.numeric(gene.universe.default)){
		stop("\tExternal Enrichment ERROR: A gene universe is required and should be numeric.\n")
	}
	if(is.na(multiple.testing.method)||!is.element(multiple.testing.method,c('bonferroni','holm','hochberg','hommel','BH','fdr','BY'))){
		stop("\tExternal Enrichment ERROR: A multiple testing method is required.\n")
	}
	if(is.na(gene.list)||!is.vector(gene.list)){
		stop("\tExternal Enrichment ERROR: A gene list as vector is requiered.\n")	
	}
	if(is.na(associations.file)||!file.exists(associations.file)){
		stop("\tExternal Enrichment ERROR: A associations file is requiered.\n")
	}
	##	
	suppressMessages(library(biomaRt))
	cat("\tGO Terms Enrichment: Downloading and formating Gene Ontology data.\n")
	GODataExternaltemp<-read.csv(file=associations.file,header=TRUE,sep="\t")
	colnames(GODataExternaltemp)<-c('GO_GO_id','GO_GO_name','GO_ensembl_Gene_ID','GO_Gene_Symbole')
	GODataExternaltemp[GODataExternaltemp == ""]<-NA
	GODataExternaltemp<-subset(GODataExternaltemp, !is.na(GODataExternaltemp$GO_ensembl_Gene_ID) || !is.na(GODataExternaltemp$GO_GO_id))
	GO.data<-GODataExternaltemp
	cat("\tGO Terms Enrichment: Gene universe: ")
	gene.universe<-gene.universe.default
	if(length(unique(GO.data$GO_ensembl_Gene_ID))>=gene.universe.default){
		gene.universe<-length(unique(GO.data$GO_ensembl_Gene_ID))	
		cat(gene.universe,"\n")
	}else{
		gene.universe<-gene.universe.default
		cat(gene.universe,"\n")
	}
	##Number of associated genes to each X in the gene universe
	total.genes.per.GO<-table(GO.data$GO_GO_id)
	total.genes.per.GO<-total.genes.per.GO[total.genes.per.GO != 0]
	##Number of associated genes to each X in the gene list
	GO.data.sub<-subset(GO.data,GO.data$GO_ensembl_Gene_ID %in% gene.list)
	list.genes.per.GO<-table(GO.data.sub$GO_GO_id)
	list.genes.per.GO<-list.genes.per.GO[list.genes.per.GO != 0]
	bigMat<-data.frame(ID=character(),NAME=character(),Fa=numeric(), Fb=numeric(), Fc=numeric(), Fd=numeric(), FPV=numeric(), GENE=character(),stringsAsFactors=FALSE)
	cat("\tGO Terms Enrichment: Statistical analysis.\n")
	for(GO in unique(GO.data.sub$GO_GO_id)){
		tempmat<-matrix(c(list.genes.per.GO[GO],(length(gene.list)-list.genes.per.GO[GO]),total.genes.per.GO[GO],(gene.universe-total.genes.per.GO[GO])),ncol = 2, nrow = 2)
		if((is.na(list.genes.per.GO[GO])==TRUE)||(is.na(length(gene.list)-list.genes.per.GO[GO]) == TRUE)||(is.na(total.genes.per.GO[GO]) == TRUE)||(is.na(gene.universe-total.genes.per.GO[GO]) == TRUE)){
			print(tempmat)
			print(GO)
			print("NA value")
			next
		}
		if((is.numeric(list.genes.per.GO[GO])==FALSE)||(is.numeric(length(gene.list)-list.genes.per.GO[GO]) == FALSE)||(is.numeric(total.genes.per.GO[GO]) == FALSE)||(is.numeric(gene.universe-total.genes.per.GO[GO]) == FALSE)){
			print(tempmat)
			print(GO)
			print("Non numeric value")
			next
		}
		if((is.infinite(list.genes.per.GO[GO])==TRUE)||(is.infinite(length(gene.list)-list.genes.per.GO[GO]) == TRUE)||(is.infinite(total.genes.per.GO[GO]) == TRUE)||(is.infinite(gene.universe-total.genes.per.GO[GO]) == TRUE)){
			print(tempmat)
			print(GO)
			print("Infinite value")
			next
		}
		if((list.genes.per.GO[GO] < 0)||((length(gene.list)-list.genes.per.GO[GO]) < 0)||(total.genes.per.GO[GO] < 0)||((gene.universe-total.genes.per.GO[GO]) < 0)){
			print(tempmat)
			print(GO)
			print("Negative value")
			next
		}
		ft<-fisher.test(tempmat,alternative="greater")
		tempmat[1,1]<-tempmat[1,1]-1
		GO.name<-unique(subset(GO.data,GO.data$GO_GO_id == GO)[,'GO_GO_name'])
		associated.genes<-paste(subset(GO.data.sub,GO.data.sub$GO_GO_id == GO)[,'GO_ensembl_Gene_ID'],subset(GO.data.sub,GO.data.sub$GO_GO_id == GO)[,'GO_Gene_Symbole'],sep=":",collapse=",")
		bigMat<-rbind(bigMat,data.frame(ID=GO,NAME=GO.name,Fa=list.genes.per.GO[GO],Fb=(length(gene.list)-list.genes.per.GO[GO]),Fc=total.genes.per.GO[GO],Fd=(gene.universe-total.genes.per.GO[GO]),FPV=ft$p.value,GENE=associated.genes))
	}
	bigMat.adjusted<-cbind(bigMat,'adjusted_pvalue'=p.adjust(as.vector(bigMat[,'FPV']),method=multiple.testing.method))
	return(bigMat.adjusted)
}
##################
##Gene Orthology##
##################
##Input: a subject bioMartID, a query bioMart ID, a gene list 
##Output: a data frame with two columns
#Comment: work only for ensembl gene ID
ortho.FUN<-function(
	subset.martID=NA,
	query.genome=NA,
	genes.list=NA,
	confidenceCT=1,
	subset.duplicates=TRUE,
	query.duplicates=TRUE
){
	cat("\tOrthology: Test on parameters.\n")
	if(is.na(subset.martID)){
		stop("\tOrthology ERROR: A subset bioMart ID is required.\n")
	}
	if(is.na(query.genome)){
		stop("\tOrthology ERROR: A query genome is required.\n")
	}
	if(!is.vector(genes.list)){
		stop("\tOrthology ERROR: A list of ensembl gene ID is required.\n")
	}
	cat("\tOrthology: Downloading and formating orthology data.\n")
	suppressMessages(library(biomaRt))
	ensembl<-useMart('ENSEMBL_MART_ENSEMBL',host="www.ensembl.org",dataset=subset.martID)
	#OR#ensembl=useMart("ensembl",dataset=martID)	#to use if the previous line rise error
	ortho.temp<-getBM(attributes=c('ensembl_gene_id',paste(query.genome,'_homolog_ensembl_gene',sep=""),paste(query.genome,'_homolog_orthology_confidence',sep="")), mart=ensembl)
	ortho<-unique(ortho.temp)
	ortho[ortho==""]<-NA
	ortho<-ortho[complete.cases(ortho),]
	colnames(ortho)<-c('subsetID','queryID','confidence')
	ortho<-subset(ortho,ortho$confidence >= confidenceCT)
	cat("\tOrthology: Formating the gene list.\n")
	replace(genes.list,genes.list=="",NA)
	genes.list<-na.omit(genes.list)
	subset.query<-subset(ortho, ortho$subset %in% genes.list)[,c('subsetID','queryID')]
	####################
	if(subset.duplicates && query.duplicates){
		cat("\tOrthology: Efficiency ",round((nrow(subset.query)*100)/length(genes.list),digits=2),"%.\n")
		return(subset.query)
	}else{
		if(! subset.duplicates){
			subset.query<-subset.query[!duplicated(subset.query$subsetID),]
		}
		if(! query.duplicates){
			subset.query<-subset.query[!duplicated(subset.query$queryID),]
		}
		cat("\tOrthology: Efficiency ",round((nrow(subset.query)*100)/length(genes.list),digits=2),"%.\n")
		return(subset.query)
	}
}

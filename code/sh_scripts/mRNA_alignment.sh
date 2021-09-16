#!/bin/bash
##mRNAseq read mapping and counting
##Dependencies: STAR, samtools, HTSeq 
##############
##Parameters
EMAIL=''
THREAD=
GENOMEDIR=''
PREFIXS=('')
DIRIN=''
GTF=''
########
##Main
#############
##Alignment##
#############
echo "Alignment"
for p in "${PREFIXS[@]}"
do
	echo $p
	OUT=$DIRIN/$p/$p.aligned.1p
	FASTQ=$(find $DIRIN/$p/ -name "*.gz" | sort | tr "\n" "," | sed 's/.$//')
	STAR --runMode alignReads --runThreadN $THREAD --genomeDir $GENOMEDIR --readFilesCommand zcat --outFileNamePrefix $OUT --outSAMtype BAM Unsorted --outSAMattributes All --outSAMmode Full --outSAMunmapped None --outFilterType Normal --outFilterMultimapNmax 1 --outSJfilterReads All --outFilterScoreMin 0 --outFilterIntronMotifs None --outSAMstrandField None --readFilesIn $FASTQ
done
echo "Second pass"
SJfiles=$(find $DIRIN -name "*SJ.out.tab" | tr "\n" "\t")
echo $SJfiles
for p in "${PREFIXS[@]}"
do
	echo $p
	OUT=$DIRIN/$p/$p.aligned.2p
	FASTQ=$(find $DIRIN/$p/ -name "*.gz" | sort | tr "\n" "," | sed 's/.$//')
	STAR --runMode alignReads --runThreadN $THREAD --genomeDir $GENOMEDIR --readFilesCommand zcat --outFileNamePrefix $OUT --outSAMtype BAM Unsorted --outSAMattributes All --outSAMmode Full --outSAMunmapped None --outFilterType Normal --outFilterMultimapNmax 1 --outSJfilterReads All --outFilterScoreMin 0 --outFilterIntronMotifs None --outSAMstrandField None --sjdbFileChrStartEnd $SJfiles --readFilesIn $FASTQ
done
wait
###############################
##Statistics on raw alignment##
###############################
echo "Statistics"
for p in "${PREFIXS[@]}"
do
	samtools flagstat $DIRIN/$p/$p.aligned.2p*.bam >> $DIRIN/$p/$p.flagstat.txt
done
for p in "${PREFIXS[@]}"
do
	samtools stats -d $DIRIN/$p/$p.aligned.2p*.bam >> $DIRIN/$p/$p.stat.txt
done
for p in "${PREFIXS[@]}"
do
	mkdir $DIRIN/$p/BAMPlot
done
for p in "${PREFIXS[@]}"
do
	plot-bamstats -p $DIRIN/$p/BAMPlot/ $DIRIN/$p/$p.stat.txt
done
#################
##Read counting##
#################
echo "Read counting"
for p in "${PREFIXS[@]}"
do
	samtools view -H $DIRIN/$p/$p.aligned.2p*.bam >> $DIRIN/$p/$p.aligned.sam&
done
for p in "${PREFIXS[@]}"
do
	samtools view $DIRIN/$p/$p.aligned.2p*.bam | cut -f18,19 --complement >> $DIRIN/$p/$p.aligned.sam
done
for p in "${PREFIXS[@]}"
do
	samtools view $DIRIN/$p/$p.aligned.sam | python htseq-count -r pos --stranded=no --mode=union - $GTF > $DIRIN/$p/$p.count
done
echo -e "mRNAseq read mapping and counting DONE\n" | cat - $DIRIN/logs.log | sendmail $EMAIL
exit

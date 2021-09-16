#!/bin/bash
##ChIPseq data read mapping and cleaning
##Depend on: bwa, samtools, GATK
##########
##Parameters
EMAIL=''
THREAD=
DIRIN=''
SAMPLE=''
GENOME=''
##########
##Main
###############
###Alignment###
###############
bwa aln -t$THREAD $GENOME $DIRIN/$SAMPLE.fastq.gz > $DIRIN/$SAMPLE.sai
bwa samse $GENOME $DIRIN/$SAMPLE.sai $DIRIN/$SAMPLE.fastq.gz | samtools view -bS - > $DIRIN/$SAMPLE.bam
###################################################
###remove unmapped reads and low mapping quality###
###################################################
samtools view -F 4 -h -q 30 -b -o $DIRIN/$SAMPLE.q1.bam $DIRIN/$SAMPLE.bam
samtools sort -O bam -o $DIRIN/$SAMPLE.q1.srt.bam $DIRIN/$SAMPLE.q1.bam
samtools index $DIRIN/$SAMPLE.q1.srt.bam
#############################
###remove duplicated reads###
#############################
java -jar MarkDuplicates INPUT=$DIRIN/$SAMPLE.q1.srt.bam OUTPUT=$DIRIN/$SAMPLE.q1.srt.rd.bam METRICS_FILE=$DIRIN/$SAMPLE.q1.srt.rd.mtr ASSUME_SORTED=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true REMOVE_SEQUENCING_DUPLICATES=true
#############################
###Remove multimapped read###
#############################
samtools view -H $DIRIN/$SAMPLE.q1.srt.rd.bam >> $DIRIN/$SAMPLE.q1.srt.rd.rm.sam
samtools view $DIRIN/$SAMPLE.q1.srt.rd.bam | grep -v -E "\sXA:" >> $DIRIN/$SAMPLE.q1.srt.rd.rm.sam
samtools view -h -b -o $DIRIN/$SAMPLE.q1.srt.rd.rm.bam $DIRIN/$SAMPLE.q1.srt.rd.rm.sam
samtools sort -O bam -o $DIRIN/$SAMPLE.q1.srt.rd.rm.srt.bam $DIRIN/$SAMPLE.q1.srt.rd.rm.bam
samtools index $DIRIN/$SAMPLE.q1.srt.rd.rm.srt.bam
######################
###Collect metrics####
######################
samtools flagstat $DIRIN/$SAMPLE.bam >> $DIRIN/$SAMPLE.bam.flagstat
samtools flagstat $DIRIN/$SAMPLE.q1.srt.bam >> $DIRIN/$SAMPLE.q1.srt.bam.flagstat
samtools flagstat $DIRIN/$SAMPLE.q1.srt.rd.bam >> $DIRIN/$SAMPLE.q1.srt.rd.bam.flagstat
samtools flagstat $DIRIN/$SAMPLE.q1.srt.rd.rm.srt.bam >> $DIRIN/$SAMPLE.q1.srt.rd.rm.srt.bam.flagstat
exit

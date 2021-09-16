#!/bin/bash
##CELseq read mapping
##Dependencies: bwa, samtools
#########
##Parameters
EMAIL=''
THREAD=
PREFIXS=('')
DIRIN=''
BWAGENOME=''
#########
##Main
for p in "${PREFIXS[@]}"
do
	echo -e "$p"
	for i in $DIRIN/$p/*r2.fastq.gz
	do
		BN=$(basename -s .fastq.gz $i)
		bwa aln -t$THREAD $BWAGENOME $i > $DIRIN/$p/$BN.sai
	done
	wait
	for i in $DIRIN/$p/*r2.sai
	do
		BN=$(basename -s .sai $i)
		bwa samse $BWAGENOME  $i $DIRIN/$p/$BN.fastq.gz | samtools view -bS - > $DIRIN/$p/$p.bam&
	done
	wait
done
wait
echo -e "CELseq read mapping DONE\n" | sendmail $EMAIL
exit

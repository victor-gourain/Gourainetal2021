#!/bin/bash
##QC of fastq file from SRA database
##Dependencies: FASTX toolkit
##########
##Parameters
EMAIL=''
DIROUT=''
DIRIN=''
########
##Main
#############
for i in $(find $DIRIN -name "*.fastq.gz");do
	#echo $i
	BN=$(basename -s .fastq.gz $i)
	echo $BN
	zcat $i | head -n 1000000 | fastx_quality_stats -Q 33 -o $DIROUT/$BN.fastx.stats
	fastx_nucleotide_distribution_graph.sh -i $DIROUT/$BN.fastx.stats -o $DIROUT/$BN.fastx.stats.nucldistr.png
	fastq_quality_boxplot_graph.sh -i $DIROUT/$BN.fastx.stats -o $DIROUT/$BN.R1.fastx.stats.qual.png
done
wait
echo -e "SRA data QC DONE\n" | sendmail $EMAIL
exit

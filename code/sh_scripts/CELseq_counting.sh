#!/bin/bash
##CELseq read counting
##Dependencies: samtools, HTSeq
#########
##Parameters
EMAIL=''
PREFIXS=('')
DIRIN=''
GTF=''
#########
##Main
for p in "${PREFIXS[@]}"
do
	samtools view $DIRIN/$p/$p.bam | python htseq-count -r pos --stranded=no --mode=union - $GTF > $DIRIN/$p/$p.count
done
wait
echo -e "CELLseq read counting DONE\n" | sendmail $EMAIL
exit

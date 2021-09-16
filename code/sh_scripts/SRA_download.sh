#!/bin/bash
EMAIL=''
DIROUT=''
SRA="$DIROUT/GSE115541_series_matrix.csv"
#############
while read r;do
	name=$(echo $r | cut -d" " -f1)
	id=$(echo $r | cut -d" " -f2)
	echo "$name : $id"
	fastq-dump --outdir $DIROUT --gzip $id
	wait
	mv $DIROUT/$id.fastq.gz $DIROUT/$name.fastq.gz
done< <(sed '1d' $SRA)
echo -e "SRA download DONE\n" | cat - $DIROUT/nohup.out | sendmail $EMAIL
exit

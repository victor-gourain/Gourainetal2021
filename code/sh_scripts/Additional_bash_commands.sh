##Peak calling with Macs2
macs2 callpeak -t <samplebamfile> -c <inputbamfile> -f BAM -B --SPMR -g mm -n <prefix> -p 1e-5 --outdir /folder/

##Trimming of adapter sequence in read with cutadapt
cutadapt -j 30 -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC -o output.fastq input.fastq.gz >> output.info

##Find transcription factor binding motif with HOMER
findMotifsGenome.pl input.bed reference_genome.fa output.folder -size <n>

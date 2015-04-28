# Script to generate biased reads
#!/bin/sh
pushd ../examples/reads_gc_bias

mkdir -p Out
reads generate null_illumina_reads_gc_bias.json
bwa mem -t 8 -p ../data/red_alga.fa.gz  Out/reads.fq > Out/temp.sam 2> /dev/null
samtools view -Sb  Out/temp.sam > Out/temp.bam 2> /dev/null
samtools sort Out/temp.bam  Out/reads
samtools index Out/reads.bam

popd
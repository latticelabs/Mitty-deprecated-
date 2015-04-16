#!/bin/bash
pushd ../examples/demo

reads --pfile=illumina_reads_low.json -v

#bwa index ../data/red_alga.fa.gz
bwa mem -t 8 -p ../data/red_alga.fa.gz  Out/null_reads_low.fq > Out/temp.sam
samtools view -Sb  Out/temp.sam > Out/temp.bam
samtools sort Out/temp.bam  Out/null_reads_low
samtools index Out/null_reads_low.bam

perfectbam --inbam=Out/null_reads_low.bam

popd
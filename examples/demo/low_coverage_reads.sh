#!/bin/bash
# This script is run by the documentation generator to create the low coverage data that is used to make the BAM analysis plots
# The 30x data comes from few chromosomes and has so many reads that it is a little slow to process

pushd ../examples/demo

reads generate illumina_reads_low.json -v

#bwa index ../data/red_alga.fa.gz
bwa mem -t 8 -p ../data/red_alga.fa.gz  Out/null_reads_low.fq > Out/temp.sam
samtools view -Sb  Out/temp.sam > Out/temp.bam
samtools sort Out/temp.bam  Out/null_reads_low
samtools index Out/null_reads_low.bam

perfectbam Out/null_reads_low.bam

popd
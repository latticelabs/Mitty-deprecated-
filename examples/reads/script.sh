#!/bin/sh
set -xe

: This shell script tests generating reads
mkdir -p Out

: Generate null reads
reads generate null_illumina_reads.json

bwa mem -t 8 -p ../data/red_alga.fa.gz  Out/reads.fq > Out/temp.sam
samtools view -Sb  Out/temp.sam > Out/temp.bam
samtools sort Out/temp.bam  Out/reads
samtools index Out/reads.bam
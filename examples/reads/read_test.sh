#!/bin/sh
set -xe

: This shell script tests generating reads
mkdir -p Out

: Generate null reads
reads --pfile=null_illumina_reads.json -v

bwa mem -t 8 -p ../data/chimera.fa.gz  Out/reads.fq > Out/temp.sam
samtools view -Sb  Out/temp.sam > Out/temp.bam
samtools sort Out/temp.bam  Out/reads
samtools index Out/reads.bam

: Generate variant database
#genomes generate --pfile=variants.json -v

: Generate reads from variant sample
#reads --pfile=var_illumina_reads.json -v

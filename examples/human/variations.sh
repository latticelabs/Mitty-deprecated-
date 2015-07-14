#!/bin/bash
: File to generate a set of SNPs
set -xe

genomes generate variations.json -v -p
reads generate illumina_reads.json -v -p
pybwa ~/Data/hg38.fa.gz Out/reads_c.fq Out/reads_c.bam -p -v

pushd `dirname $0`

PROGDIR=../../mitty
OUTDIR="Out/"
mkdir -p ${OUTDIR}

: Generate mutations
python ${PROGDIR}/denovo.py  --pfile=variations.json -v

: Show us the files
#cat ${OUTDIR}/test.vcf
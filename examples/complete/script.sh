#!/bin/bash
set -x

PROGDIR=../../mitty/
GENOMEDIR=../data
OUTDIR=Out/

mkdir -p ${OUTDIR}

python ${PROGDIR}/denovo.py  --pfile=snp.json -v
python ${PROGDIR}/vcf2reads.py  --pfile=illumina_reads.json -v

bwa mem -t 8 -p ${GENOMEDIR}/chimera.fa.gz  ${OUTDIR}/reads.fq > ${OUTDIR}/test.sam
samtools view -Sb  ${OUTDIR}/test.sam > ${OUTDIR}/temp.bam
samtools sort ${OUTDIR}/temp.bam  ${OUTDIR}/test
samtools index ${OUTDIR}/test.bam

python ${PROGDIR}/checkbam.py  --inbam ${OUTDIR}/test.bam --out_prefix ${OUTDIR}/misalign -v
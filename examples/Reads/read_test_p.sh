#!/bin/sh
set -xe

#pushd `dirname $0`

: This shell script tests generating reads
PROGDIR=../../mitty/
GENOMEDIR=../data/
OUTDIR=Out/
PFILE=seq_reads_p.json

mkdir -p ${OUTDIR}

: Create null reads from test genome
python ${PROGDIR}vcf2reads.py --pfile=${PFILE}  -V

: Do the fake alignment to check
python ${PROGDIR}reads2bam.py -p --fa_dir=${GENOMEDIR} --fastq ${OUTDIR}reads_p.fq --bam=${OUTDIR}reads_p_aligned.bam -v
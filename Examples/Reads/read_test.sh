#!/bin/sh
set -xe

pushd `dirname $0`

: This shell script tests generating reads
PROGDIR=../../mitty/
GENOMEDIR=../data/
OUTDIR="Out/"
PFILE=seq_reads.json

mkdir -p ${OUTDIR}

: Create null reads from test genome
python ${PROGDIR}vcf2reads.py --fa_dir=${GENOMEDIR}  --out=${OUTDIR}null  --pfile=${PFILE}  --block_len=1000 --master_seed=1 -V

: Do the fake alignment to check
python ${PROGDIR}reads2bam.py  --fa_dir=${GENOMEDIR} --fastq ${OUTDIR}null.fq --bam=${OUTDIR}test.bam -v
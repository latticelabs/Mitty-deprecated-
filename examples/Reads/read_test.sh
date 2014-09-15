#!/bin/sh
set -xe

#pushd `dirname $0`

: This shell script tests generating reads
PROGDIR=../../mitty/
GENOMEDIR=../data
OUTDIR=Out/

if [ $1 == "nullseq" ]
then
  PFILE=seq_reads.json
elif [ $1 == "nullillumina" ]
then
  PFILE=null_illumina_reads.json
elif  [ $1 == "varillumina" ]
then
  PFILE=var_illumina_reads_p.json
fi

mkdir -p ${OUTDIR}

: Create reads
python ${PROGDIR}vcf2reads.py --pfile=${PFILE}  -V

: Do the fake alignment to check
python ${PROGDIR}reads2bam.py -p --fa_dir=${GENOMEDIR} --fastq ${OUTDIR}reads.fq --bam=${OUTDIR}reads_aligned.bam -v
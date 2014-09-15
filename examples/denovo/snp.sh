#!/bin/bash
: File to generate a set of SNPs
set -xe

pushd `dirname $0`

PROGDIR=../../mitty
OUTDIR="Out/"
mkdir -p ${OUTDIR}

: Generate mutations
python ${PROGDIR}/denovo.py  --pfile=snp.json -v

: Show us the files
cat ${OUTDIR}/test.vcf
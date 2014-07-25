#!/bin/bash
set -xe

pushd `dirname $0`

PROGDIR=../../mitty
OUTDIR=$(mktemp -d -t WG)

: Create whole genome file from fasta files
python ${PROGDIR}/fasta2wg.py --index=wg_chimera.json --wg=${OUTDIR}/chimera.h5 --fa=${OUTDIR}/chimera.fa.gz -v
: Use the describe command to reveal the contents of the genome file
python ${PROGDIR}/fasta2wg.py describe --wg=${OUTDIR}/chimera.h5

popd
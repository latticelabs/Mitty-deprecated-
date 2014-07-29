#!/bin/bash
: File to generate a set of simple insertions
set -xe

pushd `dirname $0`

PROGDIR=../../mitty
OUTDIR=$(mktemp -d -t WG)

: Create whole genome file from fasta files
python ${PROGDIR}/fasta2wg.py --index=../fasta2wg/wg_chimera.json --wg=${OUTDIR}/chimera.h5 --fa=${OUTDIR}/chimera.fa.gz -v

: Generate mutations
python ${PROGDIR}/mutate.py  --wg=${OUTDIR}/chimera.h5 --vcf=${OUTDIR}/chimera_simple_insert_var.vcf.gz --paramfile=simple_insert.json -v

: Show us the files
cat ${OUTDIR}/chimera_simple_insert_var.vcf

#: Sort, compress and index the VCF file
#vcf-sort ${DATADIR}/chimera_simple_insert_var.vcf > ${DATADIR}/chimera_simple_insert_var_srt.vcf
#bgzip -c ${DATADIR}/chimera_simple_insert_var_srt.vcf > ${DATADIR}/chimera_simple_insert_var_srt.vcf.gz
#tabix ${DATADIR}/chimera_simple_insert_var_srt.vcf.gz

#: Generate the mutated genome
#python vcf2seq.py --ref ${DATADIR}/chimera.h5 --vcf ${DATADIR}/chimera_simple_insert_var_srt.vcf.gz --var ${DATADIR}/chimera_simple_insert.h5 -v
#
#: Use the describe command to reveal the contents of the genome file
#python fasta2wg.py describe --wg=${DATADIR}/chimera_simple_insert.h5

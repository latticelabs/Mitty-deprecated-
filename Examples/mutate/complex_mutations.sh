#!/bin/bash
: Showing full range of mutations possible using Mitty

set -xe
DIR=Examples/Mutations
DATADIR=${DIR}/Out

: Create whole genome file from fasta files
python fasta2wg.py --index=${DIR}/wg_chimera.json --wg=${DATADIR}/chimera.h5 --fa=${DATADIR}/chimera.fa.gz -v
: Use the describe command to reveal the contents of the genome file
python fasta2wg.py describe --wg=${DATADIR}/chimera.h5

: Generate mutations
python mutate.py  --wg=${DATADIR}/chimera.h5 --vcf=${DATADIR}/chimera_complex.vcf --paramfile ${DIR}/complex_mutations.json -v

: Sort, compress and index the VCF file
vcf-sort ${DATADIR}/chimera_complex.vcf > ${DATADIR}/chimera_complex_srt.vcf
bgzip -c ${DATADIR}/chimera_complex_srt.vcf > ${DATADIR}/chimera_complex_srt.vcf.gz
tabix ${DATADIR}/chimera_complex_srt.vcf.gz

: Generate the mutated genome
python vcf2seq.py --ref ${DATADIR}/chimera.h5 --vcf ${DATADIR}/chimera_complex_srt.vcf.gz --var ${DATADIR}/chimera_complex.h5 -v

: Use the describe command to reveal the contents of the genome file
python fasta2wg.py describe --wg=${DATADIR}/chimera_complex.h5

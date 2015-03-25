#!/bin/bash
: File to generate a set of SNPs
set -xe

pushd `dirname $0`

OUTDIR="Out/"
mkdir -p ${OUTDIR}

: Show us the parameters
genomes dryrun --pfile=chimera_snp.json

: Generate mutations
genomes generate --pfile=chimera_snp.json -v

: Inspect the file
genomes inspect --dbfile=Out/chimera.db

: Write out a VCF from the file
genomes write vcf 0 --dbfile=Out/chimera.db

: Write out two VCFs from the file
genomes write vcf 1 2 --dbfile=Out/chimera.db

: Write out all VCFs from the file '(How large is your disk?)'
genomes write vcf --dbfile=Out/chimera.db

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


#!/bin/sh
set -xe

: This shell script tests generating reads
mkdir -p Out

: Generate null reads
reads --pfile=null_illumina_reads.json -v

: Generate variant database
genomes generate --pfile=variants.json -v

: Generate reads from variant sample
reads --pfile=var_illumina_reads.json -v
#!/bin/bash
set -x

OUTDIR="Out/"
mkdir -p ${OUTDIR}

: Show us the parameters
# genomes dryrun --pfile=snp.json

: Generate mutations
genomes generate --pfile=snp.json -v

#: Inspect the file
# genomes inspect --dbfile=Out/chimera_no_sfs.db

#: Write out a VCF from the file
genomes write vcf 0 --dbfile=Out/chimera_no_sfs.db
#
#: Write out two VCFs from the file
#genomes write vcf 1 2 --dbfile=Out/chimera.db
#
#: Write out all VCFs from the file '(How large is your disk?)'
#genomes write vcf --dbfile=Out/chimera.db

: Take reads
reads --pfile=illumina_reads.json -v

bwa mem -t 8 -p ../data/chimera.fa.gz  ${OUTDIR}/reads.fq > ${OUTDIR}/test.sam
samtools view -Sb  ${OUTDIR}/test.sam > ${OUTDIR}/temp.bam
samtools sort ${OUTDIR}/temp.bam  ${OUTDIR}/test
samtools index ${OUTDIR}/test.bam

samtools tview ${OUTDIR}/test.bam ../data/chimera.fa.gz
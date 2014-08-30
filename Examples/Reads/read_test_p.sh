#!/bin/sh
set -xe

pushd `dirname $0`

: This shell script tests generating reads
PROGDIR=../../mitty/
GENOMEDIR=../data/
OUTDIR="Out/"
PFILE=seq_reads_p.json

mkdir -p ${OUTDIR}

: Create null reads from test genome
python ${PROGDIR}vcf2reads.py --fa_dir=${GENOMEDIR}  --out=${OUTDIR}null  --pfile=${PFILE}  --block_len=1000 --master_seed=1 -V

: Do the fake alignment to check
python ${PROGDIR}reads2bam.py -p --fa_dir=${GENOMEDIR} --fastq ${OUTDIR}null.fq --bam=${OUTDIR}test.bam -v

#bwa mem Examples/Reads/Data/haploid_ref.fa.gz Examples/Reads/Data/null_reads.fastq > Examples/Reads/Data/aligned_null.sam
#
#: Generate a VCF file with some insertions
#python mutate.py --paramfile=Examples/Reads/simple_inserts.json -v
#vcf-sort Examples/Reads/Data/simple_inserts.vcf > Examples/Reads/Data/simple_inserts_sorted.vcf
#bgzip -c Examples/Reads/Data/simple_inserts_sorted.vcf > Examples/Reads/Data/simple_inserts_sorted.vcf.gz
#tabix Examples/Reads/Data/simple_inserts_sorted.vcf.gz
#
#: Generate the variant genome
#python vcf2seq.py --ref=Examples/Reads/Data/haploid_ref.h5 --vcf=Examples/Reads/Data/simple_inserts_sorted.vcf.gz --var=Examples/Reads/Data/genome1.h5
#: Create reads from it
#python reads.py --paramfile=Examples/Reads/genome1_reads.json -v --fastq
#
#: Do some alignment to check
#bwa mem -p Examples/Reads/Data/haploid_ref.fa.gz Examples/Reads/Data/genome1_reads.fastq > Examples/Reads/Data/genome1_aligned.sam
#
#
#: Generate a VCF file with a mixture of variants
#python mutate.py --paramfile=Examples/Reads/mixed_variants.json -v
#vcf-sort Examples/Reads/Data/mixed_variants.vcf > Examples/Reads/Data/mixed_variants_sorted.vcf
#bgzip -c Examples/Reads/Data/mixed_variants_sorted.vcf > Examples/Reads/Data/mixed_variants_sorted.vcf.gz
#tabix Examples/Reads/Data/mixed_variants_sorted.vcf.gz
#
#: Generate the variant genome
#python vcf2seq.py --ref=Examples/Reads/Data/haploid_ref.h5 --vcf=Examples/Reads/Data/mixed_variants_sorted.vcf.gz --var=Examples/Reads/Data/genome2.h5
#: Create reads from it
#python reads.py --paramfile=Examples/Reads/genome2_reads.json -v --fastq
#
#: Do some alignment to check
#bwa mem -p Examples/Reads/Data/haploid_ref.fa.gz Examples/Reads/Data/genome2_reads.fastq > Examples/Reads/Data/genome2_aligned.sam

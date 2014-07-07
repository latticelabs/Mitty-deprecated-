#!/bin/sh
set -x
: This shell script tests functions leading up to reads
mkdir -p Examples/Reads/Data


: Create the whole genome file
python fasta2wg.py --index=Examples/Reads/wg_test.json --wg=Examples/Reads/Data/haploid_ref.h5 --fa=Examples/Reads/Data/haploid_ref.fa.gz -v
bwa index Examples/Reads/Data/haploid_ref.fa.gz  # Indexing needed becasue we use bwa for alignment

: Create null reads from it
python reads.py --paramfile=Examples/Reads/null_reads.json -v --fastq

: Do some alignment to check
bwa mem Examples/Reads/Data/haploid_ref.fa.gz Examples/Reads/Data/null_reads.fastq > Examples/Reads/Data/aligned_null.sam

: Generate a VCF file with some insertions
python mutate.py --paramfile=Examples/Reads/simple_inserts.json -v
vcf-sort Examples/Reads/Data/simple_inserts.vcf > Examples/Reads/Data/simple_inserts_sorted.vcf
bgzip -c Examples/Reads/Data/simple_inserts_sorted.vcf > Examples/Reads/Data/simple_inserts_sorted.vcf.gz
tabix Examples/Reads/Data/simple_inserts_sorted.vcf.gz

: Generate the variant genome
python vcf2seq.py --ref=Examples/Reads/Data/haploid_ref.h5 --vcf=Examples/Reads/Data/simple_inserts_sorted.vcf.gz --var=Examples/Reads/Data/genome1.h5
: Create reads from it
python reads.py --paramfile=Examples/Reads/genome1_reads.json -v --fastq

: Do some alignment to check
bwa mem -p Examples/Reads/Data/haploid_ref.fa.gz Examples/Reads/Data/genome1_reads.fastq > Examples/Reads/Data/genome1_aligned.sam


: Generate a VCF file with a mixture of variants
python mutate.py --paramfile=Examples/Reads/mixed_variants.json -v
vcf-sort Examples/Reads/Data/mixed_variants.vcf > Examples/Reads/Data/mixed_variants_sorted.vcf
bgzip -c Examples/Reads/Data/mixed_variants_sorted.vcf > Examples/Reads/Data/mixed_variants_sorted.vcf.gz
tabix Examples/Reads/Data/mixed_variants_sorted.vcf.gz

: Generate the variant genome
python vcf2seq.py --ref=Examples/Reads/Data/haploid_ref.h5 --vcf=Examples/Reads/Data/mixed_variants_sorted.vcf.gz --var=Examples/Reads/Data/genome2.h5
: Create reads from it
python reads.py --paramfile=Examples/Reads/genome2_reads.json -v --fastq

: Do some alignment to check
bwa mem -p Examples/Reads/Data/haploid_ref.fa.gz Examples/Reads/Data/genome2_reads.fastq > Examples/Reads/Data/genome2_aligned.sam

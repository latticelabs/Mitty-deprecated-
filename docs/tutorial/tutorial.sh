#!/bin/bash
set -ex

genomes --help
genomes show model-list

# Generate genome database
genomes generate variations.json -v -p
genomes genome-file indel --sample-name g0_s0 --max-indel 20 reddus_genomes.h5 1


# Write out one of the samples to a VCF file and view the file
genomes genome-file write-vcf --sample-name g0_s0 reddus_genomes.h5 g0_s0.vcf
head -n 20 g0_s0.vcf

# generate reads
reads generate reads.json -v -p
# awk '/([0-9]+)D/ {print}' reads_c.fq

# Align them using BWA and view the alignments
#bwa index reddus_pentalgus.fa.gz
pybwa reddus_pentalgus.fa.gz reads_c.fq bwa.bam -p -v
samtools tview -d T -p "NC_010142.1:4400" bwa.bam reddus_pentalgus.fa.gz
samtools tview -d T -p "NC_010142.1:8900" bwa.bam reddus_pentalgus.fa.gz

# Analyse the alignments
perfectbam --perfect-bam -v -v -p bwa.bam
alindel --sample-name g0_s0 --indel-range 20 bwa_per.bam reddus_genomes.h5 bwa_indel.json
cat bwa_indel.json  # Note that indel counts are for each chromosome *copy* separately
alindel_plot -f bwa_indel.json --indel-range 20 -o bwa_indel.png


# Show aligner benchmarking: do another run with a 'crappier aligner'
pybwa reddus_pentalgus.fa.gz reads_c.fq bwa_poor.bam -p -v -P -S  # No mate rescue or pairing

# Analyse the alignments
perfectbam --perfect-bam -v -v -p bwa_poor.bam
alindel --sample-name g0_s0 --indel-range 20 bwa_poor_per.bam reddus_genomes.h5 bwa_poor_indel.json
#cat bwa_poor_indel.json  # Note that indel counts are for each chromosome *copy* separately

alindel_plot -f bwa_indel.json -l 'BWA' -f bwa_poor_indel.json -l 'BWA/poor' --indel-range 20 -o combined_indel.png

alindel --sample-name g0_s1 --indel-range 20 bwa_per.bam reddus_genomes.h5 bwa_scrambled_indel.json
alindel_plot -f bwa_indel.json -l 'BWA' -f bwa_poor_indel.json -l 'BWA/poor' -f bwa_scrambled_indel.json -l 'BWA/scambled' --indel-range 20 -o combined_indel.png


# For reads solely from neighborhoods of variants
reads generate reads_variants_only.json -v -p
pybwa reddus_pentalgus.fa.gz reads_variants_only.fq reads_vo.bam -p -v

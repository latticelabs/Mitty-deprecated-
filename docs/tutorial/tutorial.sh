#!/bin/bash
set +x

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
pybwa reddus_pentalgus.fa.gz reads_c.fq reads_c.bam -p -v
samtools tview -d T -p "NC_010142.1:4400" reads_c.bam reddus_pentalgus.fa.gz
samtools tview -d T -p "NC_010142.1:8900" reads_c.bam reddus_pentalgus.fa.gz

# TODO: Show aligner benchmarking


perfectbam --perfect-bam -v -v -p reads_c.bam
alindel --sample-name g0_s0 --indel-range 20 reads_c_per.bam reddus_genomes.h5 reddus_indel.json
head -n 20 reddus_indel.json  # Note that indel counts are for each chromosome *copy* separately
alindel_plot -f reddus_indel.json --indel-range 20
exit 0


pybwa red_alga.fa.gz reads_c.fq reads_c_2.bam -p -v -P -S  # No mate rescue or pairing

perfectbam reads_c_2.bam --catreads catreads_2.h5 -v --window 10 -p
alindel catreads_2.h5 indel_2.pkl red_alga_genomes.h5 g0_s0 --indel-range 50

plot_indel -f indel.pkl -f indel_2.pkl -l 'BWA' -l 'BWA(bad)' --indel-range 50 --win 5

#genomes write vcf Out/red_alga_genomes.h5 Out/alga  --sample_name g0_s0
#alindel Out/catreads.h5 Out/indel2.pkl --vcf Out/alga_g0_s0.vcf.gz --indel 50


# For reads solely from neighborhoods of variants
reads generate reads_variants_only.json -v -p
pybwa reddus_pentalgus.fa.gz reads_variants_only.fq reads_vo.bam -p -v

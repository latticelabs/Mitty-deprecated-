#!/bin/bash
genomes generate variations.json -v -p
reads generate reads.json -v -p

#bwa index reddus_pentalgus.fa.gz
pybwa reddus_pentalgus.fa.gz reads_c.fq reads_c.bam -p -v

perfectbam reads_c.bam --catreads catreads.h5 -v --window 10 -p
alindel catreads.h5 indel.pkl red_alga_genomes.h5 g0_s0 --indel-range 50

pybwa red_alga.fa.gz reads_c.fq reads_c_2.bam -p -v -P -S  # No mate rescue or pairing

perfectbam reads_c_2.bam --catreads catreads_2.h5 -v --window 10 -p
alindel catreads_2.h5 indel_2.pkl red_alga_genomes.h5 g0_s0 --indel-range 50

plot_indel -f indel.pkl -f indel_2.pkl -l 'BWA' -l 'BWA(bad)' --indel-range 50 --win 5

#genomes write vcf Out/red_alga_genomes.h5 Out/alga  --sample_name g0_s0
#alindel Out/catreads.h5 Out/indel2.pkl --vcf Out/alga_g0_s0.vcf.gz --indel 50
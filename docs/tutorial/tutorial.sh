#!/bin/bash
genomes generate variations.json -v -p
reads generate reads.json -v -p

bwa index red_alga.fa.gz
pybwa red_alga.fa.gz reads_c.fq reads_c.bam -p -v

perfectbam reads_c.bam --catreads catreads.h5 -v --window 10 -p
alindel catreads.h5 indel.pkl --vdb red_alga_genomes.h5 --sample_name g0_s0 --indel 50

#genomes write vcf Out/red_alga_genomes.h5 Out/alga  --sample_name g0_s0
#alindel Out/catreads.h5 Out/indel2.pkl --vcf Out/alga_g0_s0.vcf.gz --indel 50
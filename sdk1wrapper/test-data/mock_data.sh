#!/bin/bash
bwa index mock.fa.gz
genomes generate -v -p variants.json
genomes genome-file write-vcf --sample-name 'g0_s0' mock_genomes.h5 mock.vcf
reads generate -v -p reads.json
pybwa mock.fa.gz mock.fq mock.bam -p -v
perfectbam mock.bam
alindel --sample-name 'g0_s0' --indel-range 20 mock_per.bam mock_genomes.h5 mock_indel.json
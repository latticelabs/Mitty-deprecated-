#!/bin/bash
: '
A script to show how mismat can be used to benchmark aligners.
Our aligners are
1. bwa mem with default settings  -> Out/reads0.bam
2. perfectbam  -> Out/reads1.bam
3. bwa mem with some settings that should knock down performance  -> Out/reads2.bam
'

pybwa ../data/red_alga.fa.gz Out/reads_c.fq Out/reads0.bam -p
perfectbam Out/reads0.bam Out/reads1.bam Out/reads0_mis.db -v -p
perfectbam Out/reads1.bam Out/reads0_p.bam Out/reads1_mis.db -v -p
pybwa ../data/red_alga.fa.gz Out/reads_c.fq Out/reads2.bam -S -P -p
perfectbam Out/reads2.bam Out/reads2_p.bam Out/reads2_mis.db -v -p

mismat Out/reads0_mis.db Out/reads1_mis.db Out/reads2_mis.db
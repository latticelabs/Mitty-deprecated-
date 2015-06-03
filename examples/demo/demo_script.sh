#!/bin/bash
# This lists most of the commands used in the quickstart guide. Help-like commands are excluded
pushd ../examples/demo

genomes generate ../examples/demo/variations.json -v -p
genomes write vcf ../examples/demo/Out/red_alga_genomes.db ../examples/demo/Out/red_alga 3


reads generate ../examples/demo/illumina_reads.json -v -p

bwa index ../data/red_alga.fa.gz
bwa mem -t 8 -p ../examples/data/red_alga.fa.gz  ../examples/demo/Out/reads_c.fq > ../examples/demo/Out/temp.sam
samtools view -Sb  ../examples/demo/Out/temp.sam > ../examples/demo/Out/temp.bam
samtools sort ../examples/demo/Out/temp.bam  ../examples/demo/Out/reads
samtools index ../examples/demo/Out/reads.bam

perfectbam ../examples/demo/Out/reads.bam ../examples/demo/Out/reads_p.bam ../examples/demo/Out/reads_mis.db -v -p
echo $'.mode line\n select * from summary limit 2;' |  sqlite3 ../examples/demo/Out/reads_mis.db

../examples/demo/low_coverage_reads.sh
plot_align circle ../examples/demo/Out/null_reads_low  --down-sample 2
plot_align matrix ../examples/demo/Out/null_reads_low

reads generate  ../examples/demo/illumina_reads_variants_only.json -v
bwa mem -t 8 -p  ../examples/data/red_alga.fa.gz  ../examples/demo/Out/vreads_c.fq > ../examples/demo/Out/temp.sam
samtools view -Sb  ../examples/demo/Out/temp.sam > ../examples/demo/Out/temp.bam
samtools sort  ../examples/demo/Out/temp.bam  ../examples/demo/Out/vreads
samtools index  ../examples/demo/Out/vreads.bam

popd
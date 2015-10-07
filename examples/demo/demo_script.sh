#!/bin/bash
# This lists most of the commands used in the tutorial guide so you can play along with the tutorial.
# Help-like commands are excluded

genomes generate variations.json -v -p
genomes write vcf ../examples/demo/Out/red_alga_genomes.db ../examples/demo/Out/red_alga 3

reads generate illumina_reads.json -v -p

bwa index ../data/red_alga.fa.gz

pushd Out
bwa mem -t 8 -p ../../data/red_alga.fa.gz  reads_c.fq > temp.sam
samtools view -Sb  temp.sam > temp.bam
samtools sort temp.bam  reads
samtools index reads.bam

perfectbam reads.bam --debugdb reads_mis.db --catreads catreads.h5 -v -p
echo $'.mode line\n select * from summary limit 2;' |  sqlite3 reads_mis.db
popd

# pybwa ../../data/red_alga.fa.gz reads.fq reads.bam -p
bucket catreads.h5 red_alga_genomes.h5 'g0_s0' s0.indel.pkl --indel=100

../examples/demo/low_coverage_reads.sh
plot_align circle ../examples/demo/Out/null_reads_low  --down-sample 2
plot_align matrix ../examples/demo/Out/null_reads_low

reads generate  ../examples/demo/illumina_reads_variants_only.json -v
bwa mem -t 8 -p  ../examples/data/red_alga.fa.gz  ../examples/demo/Out/vreads_c.fq > ../examples/demo/Out/temp.sam
samtools view -Sb  ../examples/demo/Out/temp.sam > ../examples/demo/Out/temp.bam
samtools sort  ../examples/demo/Out/temp.bam  ../examples/demo/Out/vreads
samtools index  ../examples/demo/Out/vreads.bam


genomes generate variations.json -v -p
reads generate illumina_reads.json -v -p
pybwa ../data/red_alga.fa.gz Out/reads_c.fq Out/reads.bam -p -v
#pybwa ../data/red_alga.fa.gz Out/reads_c.fq Out/reads.bam -p -v -P -S
#perfectbam Out/reads.bam --catreads Out/catreads.h5 -v -p
perfectbam Out/reads.bam --catreads Out/catreads.h5 -v --window 50 -p
alindel Out/catreads.h5 Out/indel.pkl --vdb Out/red_alga_genomes.h5 --sample_name g0_s0 --indel 50

genomes write vcf Out/red_alga_genomes.h5 Out/alga  --sample_name g0_s0
alindel Out/catreads.h5 Out/indel2.pkl --vcf Out/alga_g0_s0.vcf.gz --indel 50
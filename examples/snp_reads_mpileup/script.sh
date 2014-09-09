#!/bin/bash
set -xe

mkdir -p out
python ../../mitty/denovo.py  --fa_dir=../data/ --param_file=snp.json   --vcf=out/snp.vcf  --master_seed=1024 -v
python ../../mitty/vcf2reads.py  --fa_dir=../data/ --pfile=illumina.json  --vcf=out/snp.vcf.gz  --out=out/snp_reads  --corrupt --master_seed=1024 -v
python ../../mitty/reads2bam.py -p --fa_dir=../data/ --fastq out/snp_reads_c.fq --bam=out/test.bam -v

samtools tview -d T -p "gi|4630864|dbj|AB026117.1|:9980" out/test.bam ../data/chimera.fa.gz

bwa mem -t 8 -p ../data/chimera.fa.gz  out/snp_reads_c.fq > out/test.sam

samtools view -Sb  out/test.sam > out/temp.bam
samtools sort out/temp.bam  out/test

samtools index out/test.bam

samtools tview -d T -p "gi|4630864|dbj|AB026117.1|:9980" out/test.bam ../data/chimera.fa.gz

python ../../mitty/checkbam.py  --inbam out/test.bam --fout out/misalign.csv -v



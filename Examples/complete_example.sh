#!/bin/sh
# A simple, complete example that takes you from the sample fasta data files (and sample parameter files) to an
# aligned BAM that you can inspect

set -x  # We like to see the commands as they execute

PROGDIR="../"
DATADIR="../Data"
CHROM=1
SEQFILE="porcine_circovirus"
MUTPARFILE="mutation_par.json"
VCFFILE="variants.vcf"
READPARFILE="read_par.json"

: Convert fasta to smalla format using converta.py
python $PROGDIR/converta.py $DATADIR/$SEQFILE.fa $DATADIR/$SEQFILE.smalla

: Generate mutations and write out a VCF file
python $PROGDIR/mutate.py --chrom=$CHROM  --ref=$DATADIR/$SEQFILE.smalla  --vcf=$VCFFILE  --paramfile=$MUTPARFILE  -v
cat $VCFFILE
cat $VCFFILE.info

: Use vcf2seq to generate mutated sequence from VCF and reference sequence
python $PROGDIR/vcf2seq.py $DATADIR/$SEQFILE.smalla ex1_mutated.smalla 1 $VCFFILE.gz

: Use reads to generate a bam file full of reads
python $PROGDIR/reads.py  --ref=ex1_mutated.smalla --paramfile=$READPARFILE --coverage=10 --out=sim_reads -c

: Use cheata to fake align the reads according to the coordinates we store in the seq id
#python ../cheata.py --inbam=sim_reads.bam --outbam=aligned.bam
python ../cheata.py --inbam=sim_reads_c.bam --outbam=aligned.bam

: Use samtools to generate a VCF file of the variants
samtools mpileup -uf $DATADIR/porcine_circovirus.fa aligned.bam | bcftools view -bvcg - > var.raw.bcf
bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > mpileup.vcf

: Now compare the original VCF file with the detected one
: Original
tail -n -8 variants.vcf
: Computed by mpileup from the alignment
tail -n -8 mpileup.vcf

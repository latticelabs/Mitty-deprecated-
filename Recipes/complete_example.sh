# A complete example that takes you from the sample fasta data files (and sample parameter files) to an
# Aligned BAM that you can inspect

set -x

DATADIR="Data"
CHROM=1
SEQFILE="porcine_circovirus"
MUTPARFILE="Params/example_mutation_parameter_file.json"
VCFFILE="variants.vcf"

READPARFILE="Params/example_reads_parameter_file.json"


: Convert fasta to smalla format using converta.py
python converta.py $DATADIR/$SEQFILE.fa $DATADIR/$SEQFILE.smalla

: Generate mutations in a VCF file
python mutate.py --chrom=$CHROM  --ref=$DATADIR/$SEQFILE.smalla  --vcf=$DATADIR/$VCFFILE  --paramfile=$MUTPARFILE  -v
cat $DATADIR/$VCFFILE
cat $DATADIR/$VCFFILE.info

pushd $DATADIR
: Use samtools to compress and index the variant file
bgzip -c $VCFFILE > $VCFFILE.gz
tabix -p vcf $VCFFILE.gz
popd


: Use vcf2seq to generate mutated sequence from VCF and reference sequence
python vcf2seq.py $DATADIR/$SEQFILE.smalla $DATADIR/mutated.smalla 1 $DATADIR/$VCFFILE.gz

: Use reads to generate a bam file full of reads
python reads.py  --paramfile=$READPARFILE


: Use cheata to fake align the reads according to the coordinates we store in the seq id
pushd Data
samtools sort corrupted_reads.bam sorted_corrupted_reads
samtools index sorted_corrupted_reads.bam
popd

python cheata.py --inbam=Data/sorted_corrupted_reads.bam  --outbam=Data/cheat_alignment.bam

: Use samtools to index this fake alignment
pushd Data
samtools sort cheat_alignment.bam sorted_cheat_alignment
samtools index sorted_cheat_alignment.bam
popd


: Use samtools to create a real alignment
pushd Data
bwa index porcine_circovirus.fa
samtools bam2fq corrupted_reads.bam > raw_reads.fq
bwa mem porcine_circovirus.fa raw_reads.fq > aligned.sam
samtools view -Sb aligned.sam > aligned.bam
samtools sort aligned.bam aligned_sorted
samtools index aligned_sorted.bam
popd

: Use samtools to generate a VCF file of the variants
pushd Data
samtools mpileup -uf porcine_circovirus.fa aligned_sorted.bam | bcftools view -bvcg - > var.raw.bcf
bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > mpileup.vcf

: Now compare the original VCF file with the detected one
: Original
tail -n -8 variants.vcf
: Computed by mpileup from the alignment
tail -n -8 mpileup.vcf
popd


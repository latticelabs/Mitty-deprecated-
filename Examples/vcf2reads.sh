#Example running script to go from VCF file to reads and then align the reads and show alignment using samtools

FASTA="porcine_circovirus.fa"

pushd ../
python reads.py --ref=Data/$FASTA  --read_profile=Params/example_read_profile.py --block_len=3000 --vcf=Data/onesnp.vcf.gz
#python reads.py --ref=Data/$FASTA  --read_profile=Params/example_read_profile.py --block_len=3000
pushd Data
#bwa index $FASTA
#samtools sort raw_reads.bam sorted_raw_reads
#samtools index sorted_raw_reads.bam
samtools bam2fq raw_reads.bam > raw_reads.fq
bwa mem -p $FASTA raw_reads.fq > aligned.sam
#bwa mem $FASTA raw_reads.fq > aligned.sam
samtools view -Sb aligned.sam > aligned.bam
samtools sort aligned.bam aligned_sorted
samtools index aligned_sorted.bam

samtools tview aligned_sorted.bam $FASTA
popd
popd
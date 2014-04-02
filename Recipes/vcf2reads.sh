#Example running script to go from VCF file to reads and then align the reads and show alignment using samtools

pushd ../
python reads.py --read_profile=Params/example_read_profile.py --block_len=3000 --vcf=Data/onesnp.vcf.gz
pushd Data
#bwa index porcine_circovirus.fa
samtools sort test.bam test_sorted
samtools index test_sorted.bam
samtools bam2fq test.bam > test.fq
#bwa mem -p porcine_circovirus.fa test.fq > aligned.sam
bwa mem porcine_circovirus.fa test.fq > aligned.sam
samtools view -Sb aligned.sam > aligned.bam
samtools sort aligned.bam aligned_sorted
samtools index aligned_sorted.bam

samtools tview aligned_sorted.bam porcine_circovirus.fa
popd
popd
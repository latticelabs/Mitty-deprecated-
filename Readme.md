Mitty is a collection of modules and scripts that enable us to generate simulated genomic data to test our algorithms.
The scripts allow us to simulate mutations on a reference sequence/genome and then simulate reads from that mutated
sequence/genome.

                    mutation
                   parameters
                       |
                       V
                    --------
                   |        |
       ref seq --->| mutate |----> VCF file
                   |        |
                    --------


                     read
                   parameters
                       |
                       V
                    --------
       ref seq --->|        |
                   |        |----> reads (FASTQ)
         VCF 1 --->|        |
         VCF 2 --->| reads  |----> ideal reads (FASTQ)
             .     |        |
             .     |        |
                    --------

       This generates simulated non-ideal reads as well as ideal reads. You can pass more than one VCF file to simulate
       reads from a heterozygous or otherwise mixed sample (e.g. from tumor samples)


The main modules are:

    mutate.py - Given a reference sequence and mutation instructions generate a mutated sequence and a VCF file
    reads.py  - Given input sequence(s) and read instructions generate simulated reads as a FASTQ file(s)

Each module is designed to run as a script. Typing `python mutate.py -h` or `python reads.py -h` will list usage and
input requirements including parameter file formats. The cookbook below should also help with typical use cases.

There are two branches in the repository:

    master - stable working code
    dev    - code could be unstable/unworking but will have the latest experimental stuff going on


Cookbook
----------------

#### Generate reads given a reference sequence and a VCF file

```
cat ref.fa | vcf-consensus file.vcf.gz > out.fa #Use VCF tools to create consensus sequence
python reads.py --ref=out.fa --read_len=50 --read_count=100  #Generate the reads
```

#### Test the read generator with ideal reads from a reference genome

```
python reads.py --ref=porcine_circovirus.fa --read_len=50 --read_count=100  #Generate the reads
bwa index porcine_circovirus.fa
bwa mem porcine_circovirus.fa simulated_reads.fastq > aligned.sam
samtools view -Sb aligned.sam > aligned.bam
samtools sort aligned.bam aligned_sorted
samtools index aligned_sorted.bam
samtools tview aligned_sorted.bam porcine_circovirus.fa
```

#### Test the SNP generator with ideal reads from a mutated genome
```
python mutate.py --ref=porcine_circovirus.fa --out=mutated --paramfile=params.py  #Generate SNPs
python reads.py --ref=mutated.fa --read_len=50 --read_count=100  #Generate the reads
bwa index porcine_circovirus.fa
bwa mem porcine_circovirus.fa simulated_reads.fastq > aligned.sam
samtools view -Sb aligned.sam > aligned.bam
samtools sort aligned.bam aligned_sorted
samtools index aligned_sorted.bam
samtools tview aligned_sorted.bam porcine_circovirus.fa   # You should see the SNPs marked out
samtools mpileup -uf porcine_circovirus.fa aligned_sorted.bam | bcftools view -bvcg - > var.raw.bcf
bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf
# You should be able to compare var.filt.vcf and mutated_variants.vcf and verify they are the same
```

Where `params.py` is the same as the example parameter file in the repository

```
"""Example parameter file for mutate program

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
snp = {
  'p': 0.01
}
```

Dev notes
=========

Design choices
--------------
### Choice to output just VCF rather than a mutated sequence
The first version of the code actually output a mutated sequence AND a VCF file. `reads` then used the mutated sequence
to generate reads. This posed a scaling problem where we would end up with many GB of mutated sequences and would be an
especial problem when we do heterozygous reads (we would need to load and store in memory multiple giant sequences). The
decision was made, therefore, to simple output VCF files - the deltas as it were.

### Parameter files for mutate
1. I chose to use parameter files because we often want to rerun experiments and it became clear early on that there would
be a lot of parameters.
1. I chose to use python for the parameter file for parsimony and flexibility
1. The parameter distribution between file and command line was based on predictions of which parameters we could
experiment with most during testing
1. At this time I do not know whether having everything on the commandline would be better for PIPITOR or if param files
are preferred for Platform integration, but either way is a short code reorganization that can be done quickly at the
time of integration.

Trivia
======
Since you were dying to know: Mitty comes from James Thurber's "The Secret Life of Walter Mitty" one of my favourite
pieces from one of my favourite authors. Though [Wikipedia][wiki] has a less favourable interpretation of what Walter Mitty
stands for I follow the interpretation found in the 2013 movie of the same name. Life is difficult and full of
insurmountable obstacles. If you do not even dream that you have surmounted these obstacles how are you going to even
start?

[wiki]: http://en.wikipedia.org/wiki/The_Secret_Life_of_Walter_Mitty
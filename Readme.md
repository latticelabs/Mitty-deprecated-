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
                   |        |----> side car file with sim params
                    --------


                     read
                   parameters
                       |
                       V
                    --------
       ref seq --->|        |
                   |        |----> reads (BAM)
       VCF1.gz --->|        |
                   | reads  |----> ideal reads (BAM)
                   |        |
                   |        |----> side car file with sim params
                    --------

       This generates simulated non-ideal reads as well as ideal reads. The VCF files are compressed by
       `bgzip` and indexed by `tabix` to allow us to load variants based on positional index.


The main modules are:

    mutate.py - Given a reference sequence and mutation instructions generate a mutated sequence and a VCF file
    reads.py  - Given input sequence(s) and read instructions generate simulated reads as a FASTQ file(s)

Each module is designed to run as a script. Typing `python mutate.py -h` or `python reads.py -h` will list usage and
input requirements including parameter file formats. The cookbook below should also help with typical use cases.

There are two branches in the repository:

    master - stable working code
    dev    - code could be unstable/unworking but will have the latest experimental stuff going on

The code requires the following non-standard modules

    BioPython   - pip install biopython --user
    PyVCF       - pip install pyvcf --user
    pysam       - pip install pysam -- user


Subdirectories
--------------
    Params      - example parameter files for mutate.py and reads.py
    Recipes     - snippets of code (shell scripts and python scripts) to do/show particular tasks. useful for devs and
                  users alike
    Data        - test data for the programs



Data
----
porcine_circovirus.fa - 702bp (http://www.ncbi.nlm.nih.gov/nuccore/AY735451.1)
adenovirus.fa   -  34094bp  (http://www.ncbi.nlm.nih.gov/nuccore/AB026117.1)



Cookbook
----------------

#### Compress a VCF file with tabix




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


    bgzip test.vcf
    tabix -p vcf test.vcf.gz


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


Algorithms
----------
### Generating reads from a ref seq and VCF(s)
We want to be mindful that genomic sequences can be quite large and we might not want to carry multiple copies of such
sequences around, especially when we are simulating reads from heterogenous sequences. For this reason we choose to
avoid generating the complete sequences or ever loading the whole ref seq into memory at once.

The conceptual way to generate reads is as follows:

1. Generate the variant sequences as a whole -> there may be more than one variant sequence depending on whether there
are multiple variants at the same locus. The number of variant sequences depends on the combinatorials of all the
multiple variants
2. Take sections from each sequence at random locations, according to a predetermined distribution

If we want to avoid physically generating all combinations of variant sequences we do the following.

1. We move along the ref-seq in blocks
2. We consider only the variants within that block
3. We generate all variant sequences for that block and then generate reads randomly within that block

Notes:
1. How to handle non-local mutations? (translocation etc.)
2. How to handle multiple variants - how to link variants in one mutation with variants in another? - will currently do
a simple combinatorial with equal weights
3. How to handle the seams between blocks?
4. The number of reads from each var seq block should be adjusted to maintain coverage


V 0.2.1 (This will leave seams)
-------------------------------

                  block 0                      block 1
    Ref seq |---------------------------|----------------------------|----------------- .....


    Var seqs
            |------------------|
            |------------------------|                 --> generate reads to give required coverage
            |------------------------------|
            |--------------|


                                        |------------------|
                                        |----------------------------------|
                                        |----------------------------|
                                        |-------------------------------------|
                                        |--------------|


Notes:
1. Depending on how annoying seams are they will be fixed in V 0.2.2


Trivia
======
Since you were dying to know: Mitty comes from James Thurber's "[The Secret Life of Walter Mitty][mitty]" one of my
favourite pieces from one of my favourite authors. Though [Wikipedia][wiki] has a less favourable interpretation of what
Walter Mitty stands for, I always thought that Thurber was celebrating the dreamer within each us who spices up the banal
parts of life with a little harmless fantasy.

[mitty]: http://www.newyorker.com/archive/1939/03/18/390318fi_fiction_thurber?currentPage=all
[wiki]: http://en.wikipedia.org/wiki/The_Secret_Life_of_Walter_Mitty
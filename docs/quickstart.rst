Quick-start I
=============

Now that we have an overview of how Mitty is structured and some of the conventions used let us see how we can use
Mitty for a data simulation problem. Let us say our colleague Heng Li has written an aligner called `bwa` and we want
to figure out how well it does at aligning.

Our plan is as follows:

* Sprinkle some variants on the reference genome to create a sample genome
* Take 30x Illumina like reads from this sample genome
* Use `bwa` to align the reads to the reference genome
* Use `samtools mpileup` to create a VCF file to see if we found all the variants

*For this exercise we won't use the human chromosome but rather some test data included with the Mitty source
distribution in the `examples` directory.*

Adding variants
---------------
The example reference genome has four "chromosomes" let's add SNPs to the first two. We will need the `denovo` program
for this. Let's look into what kind of plugins we have available:

.. command-output::  python ../mitty/denovo.py plugins

Ok, the `snp` plugin sounds promising, what kind of parameters does it need?

.. command-output::  python ../mitty/denovo.py explain snp

Ok, we also want a refresher on how to write the parameter file (and what the command line options are) for denovo:

.. command-output::  python ../mitty/denovo.py

From this we first create a parameter file, let's call it `snp.json`:

.. literalinclude:: ../examples/snp_reads_mpileup/snp.json
    :language: json

Let's run this command and create a new genome (This complete example can be found under `/examples/snp_reads_mpileup`).
(*We only use the -v option to see what's going on*)

.. command-output:: mkdir -p ../examples/snp_reads_mpileup/out
.. command-output::  python ../mitty/denovo.py  --fa_dir=../examples/data/ --param_file=../examples/snp_reads_mpileup/snp.json   --vcf=../examples/snp_reads_mpileup/out/snp.vcf  --master_seed=1024 -v

Let's take a peek at the produced vcf

.. command-output:: cat ../examples/snp_reads_mpileup/out/snp.vcf

Things seem to have run satisfactorily, now lets generate a bag of reads from this genome.

Taking reads
------------
Let's first figure out what kind of read models we have available to us.

.. command-output::  python ../mitty/vcf2reads.py plugins

`simple_illumina` sounds promising, what kind of parameter file does it need?

.. command-output::  python ../mitty/vcf2reads.py explain simple_illumina

What parameters do we need to call `vcf2reads` and how do we structure the parameters file?

.. command-output:: python ../mitty/vcf2reads.py

From this we first create a parameter file, let's call it `illumina.json`:

.. literalinclude:: ../examples/snp_reads_mpileup/illumina.json
    :language: json

Let's run this command and create some reads from the variant genome

.. command-output::  python ../mitty/vcf2reads.py  --fa_dir=../examples/data/ --pfile=../examples/snp_reads_mpileup/illumina.json  --vcf=../examples/snp_reads_mpileup/out/snp.vcf.gz  --out=../examples/snp_reads_mpileup/out/snp_reads  --corrupt --master_seed=1024 -v


We can use `reads2bam.py` to create a cheat alignment ...

.. command-output::  python ../mitty/reads2bam.py -p --fa_dir=../examples/data/ --fastq ../examples/snp_reads_mpileup/out/snp_reads_c.fq --bam=../examples/snp_reads_mpileup/out/test.bam -v


... and view it using `samtools tview` (We, of course, need `samtools` installed for this to work)

.. command-output::  samtools tview -d T -p "gi|4630864|dbj|AB026117.1|:9980" ../examples/snp_reads_mpileup/out/test.bam ../examples/data/chimera.fa.gz

We can see one of the SNPs we put in.


TODO: Creating an alignment
---------------------
We can now use bwa to generate the alignment instead (make sure you have generated the index for the reference):

.. command-output::  bwa mem -t 8 -p ../examples/data/chimera.fa.gz  ../examples/snp_reads_mpileup/out/snp_reads_c.fq > ../examples/snp_reads_mpileup/out/test.sam
.. We pipe stderr to dev/null because we don't want it confusing Sphinx
.. sdfs  command-output::  samtools view -Sb  ../examples/snp_reads_mpileup/out/test.sam > ../examples/snp_reads_mpileup/out/test.bam
.. Need to work out why this does not work

And again check the alignment

.. command-output::  samtools tview -d T -p "gi|4630864|dbj|AB026117.1|:9980" ../examples/snp_reads_mpileup/out/test.bam ../examples/data/chimera.fa.gz

TODO: Add section of validating alignment - incoporate tool formally


Quick-start II
==============

We figured out how to

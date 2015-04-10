Quick-start I
=============

(For an overview of how Mitty is structured and some of the conventions used please see the detailed documentation.
This quickstart is for the impatient - like me - who like to type in commands and see stuff happen and go from
there)

We will introduce Mitty via a data simulation problem. Let us say our colleague Heng Li has written an aligner
called ``bwa`` and a variant caller called ``samtools`` and we want to figure out how well they do using synthetic data.

Our plan is as follows:

* Consider a small reference genome (we will use the Red Alga, *Cyanidioschyzon merolae strain 10D*, record_, ftp_)
* Sprinkle some variants on the reference genome to create a population of 1000 sample genomes (why not?)
* Take 30x Illumina like reads from one of the sample genomes
* Use `bwa` to align the reads back to the reference genome
* Assess how well `bwa` aligned the reads


.. _record: http://www.ncbi.nlm.nih.gov/assembly/GCF_000091205.1/#/def
.. _ftp: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000091205.1_ASM9120v1


Creating genomes
----------------
The example reference genome has 20 chromosomes. We'd like to add SNPs, insertions and deletions to chromosome 1 and 2.
We will need the ``genomes`` program for this

.. command-output::  genomes

We need the ``genomes generate`` tool, and we need a refresher on how to write a simulation parameter file

.. command-output:: genomes explain parameters

What are the variant models available to us?

.. command-output:: genomes list variantmodels

What kind of parameters do they need?

.. command-output::  genomes explain variantmodel snp

What are the population models available to us?

.. command-output:: genomes list populationmodels

An example parameter snippet?

.. command-output:: genomes explain populationmodel double_exp

Let us create a parameter file for the simulations, calling it ``variations.json``:

.. literalinclude:: ../examples/demo/variations.json
    :language: json

Before we run the actual simulation, let's make sure the site frequency spectrum looks satisfactory

.. command-output:: genomes dryrun --pfile ../examples/demo/variations.json


Let's run this command and create a database of simulated genomes

.. command-output:: mkdir -p ../examples/demo/Out
.. command-output:: genomes generate --pfile=../examples/demo/variations.json

(*Using the -v option will give detailed logger messages and a progress bar*).

Let's take a peek at the produced genomes database using the inspect function

.. command-output:: genomes inspect --dbfile=../examples/demo/Out/red_alga_genomes.db

We can also take a look at the site frequency spectrum in the generated population, for chromosome 1, for example

.. command-output:: genomes inspect sfs 1 --dbfile=../examples/demo/Out/red_alga_genomes.db

Right now, all our genomes are in a database. The rest of the world works in VCF files, so let's write out one of the
samples as a VCF file

.. command-output:: genomes write vcf --dbfile=../examples/demo/Out/red_alga_genomes.db 3
.. command-output:: head -n 20 ../examples/demo/Out/red_alga_genomes.db_s000003.vcf

Things seem to have run satisfactorily. Note that we produce phased VCF files and we also use the notation `1|0` since we
have complete knowledge of the phasing. Now let's generate a bag of reads from this genome.

Taking reads
------------
We will use the `reads` program to generate a fastq file of reads from the selected genome.

.. command-output:: reads

We know how to get help for writing the parameter file

.. command-output:: reads explain parameters

What sort of read models are currently available?

.. command-output:: reads list

``simple_illumina`` sounds promising, what kind of parameter file does it need?

.. command-output:: reads explain model simple_illumina

Ok, let's put together a parameter file for our experiment, let's call it `illumina_reads.json`:

.. literalinclude:: ../examples/demo/illumina_reads.json
    :language: json

Let's run this command and create some reads from the variant genome

.. command-output::  reads --pfile=../examples/complete/illumina_reads.json


Testing alignment accuracy of BWA MEM
-------------------------------------
First, use BWA-MEM to create an alignment:

.. command-output:: bwa index ../examples/data/red_alga.fa.gz
.. command-output:: bwa mem -t 8 -p ../examples/data/red_alga.fa.gz  ../examples/demo/Out/vreads_c.fq > ../examples/demo/Out/temp.sam
    :shell:
.. command-output::  samtools view -Sb  ../examples/demo/Out/temp.sam > ../examples/demo/Out/temp.bam
    :shell:

.. command-output:: samtools sort ../examples/demo/Out/temp.bam  ../examples/demo/Out/vreads
.. command-output:: samtools index ../examples/demo/Out/vreads.bam

We can use an alignment browser, such as Tablet, to see the alignments

.. image:: _static/vreads.bam_tablet.png

Note how the reads are restricted to the variant locations, as we asked for.


Then use the benchmarking tool `perfectbam` to analyze alignment performance and create a perfectly aligned BAM file

.. command-output:: perfectbam --inbam=../examples/demo/Out/vreads.bam



Creating a cheat alignment
--------------------------
We can use `reads2bam.py` to create a cheat alignment ...

.. command-output::  reads2bam -p --fa_dir=../examples/data/ --fastq ../examples/complete/out/reads.fq --bam=../examples/complete/out/reads.bam -v


... and view it using `samtools tview` (We, of course, need `samtools` installed for this to work)

.. command-output::  samtools tview -d T -p "gi|4630864|dbj|AB026117.1|:6950" ../examples/complete/out/reads.bam ../examples/data/chimera.fa.gz

We can see one of the SNPs we put in!


Creating an alignment
---------------------
We can now use bwa to align the simulated reads (make sure you have generated the index for the reference):

.. command-output::  bwa mem -t 8 -p ../examples/data/chimera.fa.gz  ../examples/complete/out/reads.fq > ../examples/complete/out/test.sam
    :shell:

.. command-output::  samtools view -Sb  ../examples/complete/out/test.sam > ../examples/complete/out/temp.bam
    :shell:

.. command-output:: samtools sort ../examples/complete/out/temp.bam  ../examples/complete/out/test
.. command-output:: samtools index ../examples/complete/out/test.bam


And the check the alignment

.. command-output::  samtools tview -d T -p "gi|4630864|dbj|AB026117.1|:6950" ../examples/complete/out/test.bam ../examples/data/chimera.fa.gz

We can see one of the SNPs we put in!


Assessing alignment accuracy
----------------------------
The ``checkbam`` tool allows us to read in a ``.bam`` file generated by an aligner and check alignment accuracy provided
the original reads were produced by Mitty.

.. command-output:: checkbam  --inbam ../examples/complete/out/test.bam --out_prefix ../examples/complete/out/ -v

This is the ``.json`` file with a summary of alignment performance

.. command-output:: cat ../examples/complete/out/misalign.json


and this is the ``.csv`` file with details of how the reads were misaligned

.. command-output:: head -n 10 ../examples/complete/out/misalign.csv


Quick-start II
==============

TODO: Population simulations tutorial
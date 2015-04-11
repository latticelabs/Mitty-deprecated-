Quick-start
===========

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
    :linenos:

Let's run this command and create some reads from the variant genome

.. command-output::  reads --pfile=../examples/complete/illumina_reads.json

Testing alignment accuracy of BWA MEM
-------------------------------------
First, use BWA-MEM to create an alignment:

.. command-output:: bwa index ../examples/data/red_alga.fa.gz
.. command-output:: bwa mem -t 8 -p ../examples/data/red_alga.fa.gz  ../examples/demo/Out/reads_c.fq > ../examples/demo/Out/temp.sam
    :shell:
.. command-output:: samtools view -Sb  ../examples/demo/Out/temp.sam > ../examples/demo/Out/temp.bam
    :shell:

.. command-output:: samtools sort ../examples/demo/Out/temp.bam  ../examples/demo/Out/reads
.. command-output:: samtools index ../examples/demo/Out/reads.bam

We can use an alignment browser, such as [Tablet]_, to see the alignments

.. image:: _static/reads_bam_tablet.png

(You can see that the coverage is what we dialled in. The one red colored base is a read error - since we are inspecting
the corrupted reads file)

.. [Tablet] Milne I, Stephen G, Bayer M, Cock PJA, Pritchard L, Cardle L, Shaw PD and Marshall D. 2013. Using Tablet for visual exploration of second-generation sequencing data. Briefings in Bioinformatics 14(2), 193-202. (`web`_)
.. _web: http://ics.hutton.ac.uk/tablet/


Analyse the alignment
.....................

We use ``perfectbam`` to both analyze the performance of the aligner and to create a perfectly aligned ``.bam`` file from the
``.bam`` file output by the aligner.

.. command-output::  perfectbam --inbam=../examples/demo/Out/reads.bam

`perfectbam` creates a `.json` file with some summary performance statistics and a `.csv` file with details about the
misaligned reads.

.. command-output:: head -n 15 ../examples/demo/Out/reads_summary.json

.. command-output:: head -n 10 ../examples/demo/Out/reads_misaligned.csv

We all know how exciting it is to scroll through a .csv file, so we have some additional tools to analyse the misaligned
reads and display the information graphically (These require Matplotlib_ to be installed).

.. command-output:: plot_align circle ../examples/demo/Out/reads.bam
.. #image:: _static/red_alga_circle_plot_whole.png
.. image:: ../examples/demo/Out/reads_circle_plot.png

This is a plot of where the reads end up after alignment. As can be seen, while a small percentage of reads end up in
different chromosomes, most of the misalignments are local (the lines end up on the same chromosome).

.. command-output:: plot_align matrix ../examples/demo/Out/reads.bam
.. image:: ../examples/demo/Out/reads_matrix_plot.png


.. _Matplotlib: http://matplotlib.org/index.html

Taking reads just from regions around variants
----------------------------------------------
By setting the ``variants_only`` flag to *True* in the parameters file we can selectively generate reads from regions bordering
a variation.

.. AArgh. Could not get Sphinx to
..  a) highlight the line (https://github.com/ClusterHQ/flocker/issues/337) - I think the lame green highlight matches the background
..  b) show matching line numbers (:lineno-match: and :lineno-start: make the code block disappear)

.. literalinclude:: ../examples/demo/illumina_reads_variants_only.json
    :language: json
    :emphasize-lines: 15
    :lineno-match:
    :lines: 13-17

.. command-output:: reads --pfile=../examples/demo/illumina_reads_variants_only.json
.. command-output:: bwa mem -t 8 -p ../examples/data/red_alga.fa.gz  ../examples/demo/Out/vreads_c.fq > ../examples/demo/Out/temp.sam
    :shell:
.. command-output::  samtools view -Sb  ../examples/demo/Out/temp.sam > ../examples/demo/Out/temp.bam
    :shell:

.. command-output:: samtools sort ../examples/demo/Out/temp.bam  ../examples/demo/Out/vreads
.. command-output:: samtools index ../examples/demo/Out/vreads.bam

We can use an alignment browser, such as Tablet, to see the alignments

.. image:: _static/vreads_bam_tablet.png

Note how the reads are restricted to the variant locations, as we asked for.


Analysing variant callers
-------------------------
TODO
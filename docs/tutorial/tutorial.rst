Tutorial
========

We will introduce Mitty via a data simulation problem. Let us say our colleague Heng Li has written an aligner
called ``bwa`` and we want to figure out how well it does using synthetic data.

Our plan is as follows:

* Consider a small reference genome (we will use a made up genome *Reddus pentalgus* made up of the five largest
  chromosomes of the Red Alga, *Cyanidioschyzon merolae strain 10D*, record_, ftp_ and pretend it's a diploid organism)
* Sprinkle some variants on the reference genome to create a population of 1000 sample genomes (why not?)
* Take 10x Illumina like reads from one of the sample genomes
* Use `bwa` to align the reads back to the reference genome
* Assess how well `bwa` aligned the reads


.. _record: http://www.ncbi.nlm.nih.gov/assembly/GCF_000091205.1/#/def
.. _ftp: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000091205.1_ASM9120v1

--------------

For the impatient, the commands we'll be using for this experiment are:

.. literalinclude:: tutorial.sh

(This file is included in the source distribution under ``docs/tutorial/tutorial.sh``)

--------------

Creating genomes
----------------
The example reference genome has 5 chromosomes. We'd like to add SNPs, insertions and deletions to all the chromosomes.
Since this is for testing purposes we would like to have a large number of indels distributed across a range of lengths,
rather than the power-law distribution we see in nature (which gives us many short indels and far fewer longer indels).

To perform this simulation, we need to create a parameter file describing our simulation requirements and pass this to
the ``genomes`` program. In the parameter file we need to describe three components:

 * the individual variant models
 * the site frequency spectrum model and,
 * the population generation model

The variant models are used to sprinkle variations (mutations) on top of an input reference genome. These variations are
put into a master list of variations. Each stock variation model has a parameter, ``p``, that sets the
per base probability of observing the variation in a sample. For example, setting ``p=0.001`` for the SNP plugin
and ``p=0.0001`` for the delete plugin will give us SNP and delete densities that are similar to that observed in typical
human samples.

Internally, the plugins scale the probability according to the site frequency spectrum such that samples
drawn from the master list contain variations obeying the site frequency spectrum and, overall, the desired ``p`` value.

The population generation model specifies how to construct individual samples from the master list of variants. The
default population model randomly selects variations for each sample from the master list following the site frequency
spectrum. This default model is sufficient for large population simulations.

For greater detail please see [theory and algorithms].

The genomes program performs several tasks. Typing ``genomes`` on the command line will produce a list of sub-commands.

.. command-output::  genomes

For the task at hand we need the ``genomes generate`` tool. We need to construct a parameter file which contains,
among other things, specifications for variant models, spectrum models and population models.

.. command-output:: genomes show parameters

What are the models available to us?

.. command-output:: genomes show model-list

What kind of parameters does, for example, the SNP model need?

.. command-output::  genomes show variant-model snp

The spectrum model?

.. command-output::  genomes show spectrum-model double_exp

The standard population model?

.. command-output:: genomes show population-model standard

Let us create a parameter file for the simulations, calling it ``variations.json``:

.. literalinclude:: variations.json
    :language: json
    :linenos:

*Note that if we are happy with a default plugin parameter, we can omit that in the parameter file.*

Before we run the actual simulation, let's make sure the site frequency spectrum looks satisfactory

.. command-output:: genomes generate --dry-run variations.json
    :cwd: .

The first histogram shows the site frequency spectrum that the population model will follow and the second table shows
the expected fraction of variants exhausted from the unique variants pool as more samples are generated.

Let us create a database of simulated genomes

.. command-output:: genomes generate variations.json -v
    :cwd: .

(*Using the -v option will give detailed logger messages, using the -p option will give a progress bar*).

Let's take a peek at the produced genome population database:

.. command-output:: genomes genome-file summary reddus_genomes.h5
    :cwd: .

We can also take a look at the site frequency spectrum in the generated population, for chromosome 1, for example

.. command-output:: genomes genome-file sfs reddus_genomes.h5 1
    :cwd: .

Right now, all our genomes are in a database. The rest of the world works in VCF files, so let's write out one of the
samples as a VCF file

.. command-output:: genomes genome-file write-vcf reddus_genomes.h5 g0_s0.vcf --sample-name 'g0_s0'
    :cwd: .

.. command-output:: head -n 20 g0_s0.vcf
    :cwd: .
    :shell:

Things seem to have run satisfactorily. Note that we produce phased VCF files and we also use the notation `1|0` since we
have complete knowledge of the phasing. Now let's generate a bag of reads from this genome.

Taking reads
------------
We will use the `reads` program to generate a fastq file of reads from the selected genome. We could do with a refresher
on how to create a parameters file for `reads`:

.. command-output:: reads show parameters

And what sort of read models are currently available?

.. command-output:: reads show model-list

``simple_illumina`` sounds promising, what kind of parameter file does it need?

.. command-output:: reads show read-model simple_illumina

Ok, let's put together a parameter file for our experiment, let's call it `reads.json`:

.. literalinclude:: reads.json
    :language: json
    :linenos:

Some things to note in this file are

 - The use of the flag `gzipped`. Setting this to `true` writes out gzipped FASTQ files at the expense of some runtime
   speed
 - The use of a range restriction on reads from chromosome 3 `[3, 0.4, 0.6]` which restricts reads from chrom 3 to be
   within 40% and 60% of the chromosome
 - setting of `variants_only` to be `false`. When set to `true` this allows us to get reads from only the variant
   regions. We will see an example later. This is useful for when we are doing testing and want to reduce the total
   number of reads we're processing.

Let's run this command and create some reads from the variant genome

.. command-output::  reads generate reads.json -v
    :cwd: .

The `qname` of the reads is used to encode information about the original location of the read and it's expected CIGAR.
This information is used by analysis tools to analyze aligner accuracy on these reads. The qname is convenient since it
can be stored in the FASTQ itself and seems to survive processing by all aligners tested

.. command-output::  reads --qname

.. command-output:: awk '/([0-9]+)D/ {print}' reads_c.fq | head -n 4
    :cwd: .
    :shell:


Benchmarking alignment accuracy
-------------------------------
First, use BWA-MEM to create an alignment :

.. program-output:: pybwa reddus_pentalgus.fa.gz reads_c.fq bwa.bam -p -v
    :cwd: .

.. command-output:: samtools tview -d T -p "NC_010142.1:4400" bwa.bam reddus_pentalgus.fa.gz
    :cwd: .

As part of a comparison we are going to do later on, we will run a crippled version of BWA that does not take into
account paired end information.

.. program-output:: pybwa reddus_pentalgus.fa.gz reads_c.fq bwa_poor.bam -p -v -P -S
    :cwd: .


Analyse the alignment
.....................

We use the ``perfectbam`` tool to analyze the alignment of the reads.

.. command-output:: perfectbam --help
.. command-output:: perfectbam --tags

.. command-output::  perfectbam --perfect-bam -v -v bwa.bam
    :cwd: .

.. command-output::  perfectbam bwa_poor.bam
    :cwd: .

Circle/matrix confusion plots
+++++++++++++++++++++++++++++

We can inspect the mis-alignments using a circle or matrix plot:

.. command-output:: misplot --help
.. command-output:: misplot bwa_bad.bam --matrix bwa_mat.png --circle bwa_cir.png
    :cwd: .

.. command-output:: misplot bwa_poor_bad.bam --matrix bwa_poor_mat.png --circle bwa_poor_cir.png
    :cwd: .

Note how for the normal BWA most of the confusion is local with fewer reads being scattered across the genome

.. image:: bwa_mat.png

as compared to the crippled BWA which, without the benefit of read pairing information does much more poorly

.. image:: bwa_poor_mat.png

Circle plots are fun to look at, but not as informative as the matrix plots

.. image:: bwa_cir.png
.. image:: bwa_poor_cir.png


Indel parametrized performance
++++++++++++++++++++++++++++++

We pass the perfect BAM produced in the previous stage to the ``alindel`` tool to analyze the alignment performance
parametrized by indel length.

.. command-output:: alindel --help
.. command-output:: alindel --sample-name g0_s0 --indel-range 20 bwa_per.bam reddus_genomes.h5 bwa_indel.json
    :cwd: .

.. command-output:: alindel --sample-name g0_s0 --indel-range 20 bwa_poor_per.bam reddus_genomes.h5 bwa_poor_indel.json
    :cwd: .

.. command-output:: head -n 20 bwa_indel.json
    :cwd: .

We all know how exciting it is to read `.json` files, so we have an additional tools to analyse the misaligned
reads and display the information graphically (These require Matplotlib_ to be installed).

.. command-output:: alindel_plot --help

.. command-output:: alindel_plot -f bwa_indel.json -l 'BWA' -f bwa_poor_indel.json -l 'BWA/poor' --indel-range 20 -o combined_indel.png
    :cwd: .

.. image:: combined_indel.png


.. _Matplotlib: http://matplotlib.org/index.html

Taking reads just from regions around variants
----------------------------------------------
By setting the ``variants_only`` flag to *True* in the parameters file we can selectively generate reads from regions bordering
a variation.

.. AArgh. Could not get Sphinx to highlight the line (https://github.com/ClusterHQ/flocker/issues/337) - I think the lame green highlight matches the background

.. literalinclude:: reads_variants_only.json
    :language: json
    :emphasize-lines: 15
    :lineno-match:
    :lines: 13-17

Analysing variant callers
-------------------------
TODO
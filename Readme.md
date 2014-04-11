Mitty is a collection of modules and scripts that enable us to generate simulated genomic data to test our algorithms.
The scripts allow us to simulate mutations on a reference sequence/genome and then simulate reads from that mutated
sequence/genome.

Quickstart
----------

The process for running Mitty components to create reads from a mutated genome starting from only a reference
sequence is illustrated schematically below. Please see `complete_example.sh` under the `Recipes` directory
(and the relevant parameter files under the `Params` directory) to understand the details of the command line
invocations and parameter files. For advanced users please see the rest of the docs and the `Plugins` directory
to understand how to write python code to simulate different kinds of mutations and reads.

                    ----------
         fasta     |          |---> smalla file
          file --->| converta |
                   |          |---> heada file
                    ----------

For efficiency purposes we strip the original fasta file of the header and all new lines. The
       resulting file is called a smalla file and is what the rest of the tools use. This conversion
       can be done easily using the converta.py script. The header is saved into a .smalla.heada file


                    mutation
                   parameters
                       |
                       V
                    --------
                   |        |
       ref seq --->| mutate |----> VCF file
                   |        |----> side car file with sim params
                    --------

Given a set of mutation instructions we can use mutate.py to generate a VCF file. For further
       processing the vcf file should be compressed using `bgzip` and indexed by `tabix`.



                    ---------
       ref seq --->|         |
                   | vcf2seq | ---> mutated seq
       VCF     --->|         |
                    ---------

Using the vcf2seq tool we can write out the mutations indicated by VCF into a complete mutated
       sequence saved as a smalla file.


                     read
                   parameters
                       |
                       V
                    --------
                   |        |----> corrupted ("real") reads (BAM/FASTQ)
                   |        |
           seq --->| reads  |----> ideal reads (BAM/FASTQ)
                   |        |
                   |        |----> side car file with sim params
                    --------

The reads tool enables us to take a sequence and generate simulated reads from it. The reads can
       simulate various error and property profiles of different sequencers





Each module is designed to run as a script. Typing `python mutate.py -h` or simply `python mutate.py` etc. will list
usage and input requirements.

There are two branches in the repository:

    master - stable working code
    dev    - code could be unstable/unworking but will have the latest experimental stuff going on

The code requires the following non-standard modules

    BioPython   - pip install biopython --user
    PyVCF       - pip install pyvcf --user
    pysam       - pip install pysam -- user

The code requires the following external tools to run the examples

    bgzip       - to compress VCF file
    tabix       - to index VCF file
    samtools    - indexing fasta and converting between bam and sam etc
    bwa         - alignments etc.



Subdirectories
--------------
    Params      - example parameter files for mutate.py and reads.py
    Recipes     - snippets of code (shell scripts and python scripts) to do/show particular tasks. useful for devs and
                  users alike
    Data        - test data for the programs
    Plugins     - directory where simulation models are stored

Data
----
 - porcine_circovirus.fa - 702bp (http://www.ncbi.nlm.nih.gov/nuccore/AY735451.1)
 - adenovirus.fa   -  34094bp  (http://www.ncbi.nlm.nih.gov/nuccore/AB026117.1)


Dev notes
=========

Plugin system for simulation models
-----------------------------------

Mitty comes with some prebuilt models for variant and read simulation. These models are simply Python modules that
expose a few key functions that `mutate.py` and `reads.py` use to determine variant and read characteristics. Mitty can
be extended by writing Python modules exposing these key functions and placing them in the Plugins directory.
Mitty knows to load the proper plugins simply from their names, which are mentioned in the parameters files required
for each simulation. For details please see the Readme file in the Plugins directory.

* You should use the rng passed by Mitty to ensure reproducability


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
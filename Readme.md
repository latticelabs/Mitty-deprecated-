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
for each simulation.

### Random number generators, blocked computation etc.
You will note from the stock plugins that they take in seeds for as many random number generators they use internally.
This is important to ensure reproducibility. The algorithms are designed such that the block size does not affect the
random number generation and multiple random number generators are used to avoid unexpected interactions between
different parameters, such as insertion length and position which may arise if the same random number generator is used
for both.


Algorithms
----------
### Variants, reads and CIGARS
One big goal of Mitty is to serve up realistic test data for bioinformatics algorithms, from aligners to variant
callers. Testing whether a variant caller is correctly working on the simulated data is relatively easy: we simply
compare the variant caller's VCF file with the answer book VCF generated by Mitty. It is however slightly more involved
to deduce if an aligner is correctly aligning the simulated reads and to diagnose how the performance of an aligner is
affecting the accuracy of a variant caller. To this end Mitty has a system to compute the correct read position and
CIGAR for each simulated read. This information is stored in the read's qname string so that it is easily accessible
to diagnostic programs.

`vcf2seq.py` takes in a VCF file and a reference sequence and generates three files

    mut_seq - the complete, literal, mutated sequence
    pos     - a uint16 array the same size as mut_seq, containing coordinate information that can be used to generate
              the POS value for a correctly aligned read
    diffpos - a uint16 array which is the diff of the pos array and that can be used to generate correct CIGAR strings
              for reads.

`reads.py` considers a `mut_seq` file and the supplied pos and diffpos files. If none are supplied `reads.py` assumes
that the input sequence is the reference sequence.

For every read generated from `mut_seq` the `reads.py` algorithm sets the POS (defined as the position of the first
sequence match) as the Pos value corresponding to the first Diffpos with value 1. The CIGAR is generated as follows:

Starting from the beginning of the read increment three counters M,D,I based on the Diffpos values encountered.

If the value is 1 increment the M counter, flush other counters
If the value is 0 increment the I counter and flush other counters
If the value is > 0 then increment the M counter, flush it and then increment the D counter by 1 less than the value

Flushing a counter involves writing out 2M or 3D or 5I into the cigar string if the counter value is non zero.

The algorithm is best described through a series of examples. In the examples the reference sequence is labelled `R` and
the mutated sequence is labeled `M`.



    (co-ord)     1   2   3   4   5   6   7   8   9  10  11  12
    R          | A | C | T | G | G | T | C | A | A | T | C | G |

#### Insertion

    (co-ord)     1   2   3   4                   5   6   7   8       9  10  11  12
    R          | A | C | T | G |               | G | T | C | A |   | A | T | C | G |
    M          | A | C | T | G | T | T | T | T | G | T | C | A | T | A | T | C | G |
    Pos          1   2   3   4   5   5   5   5   5   6   7   8   9   9  10  11  12
    Diffpos      1   1   1   1   0   0   0   0   1   1   1   1   0   1   1   1   1
    read 1     | A | C | T | G |                                                     POS = 1 CIGAR = '4M'
    read 2         | C | T | G | T |                                                 POS = 2 CIGAR = '3M1I'
    read 3                     | T | T | T | T |                                     POS = 0 (Unmapped read)
    read 4                         | T | T | T | G |                                 POS = 5 CIGAR = '3I1M'
    read 5                                                 | A | T | A | T |         POS = 8 CIGAR = '1M1I2M'


#### Deletion

    (co-ord)     1   2   3   4   5   6   7   8   9  10  11  12
    R          | A | C | T | G | G | T | C | A | A | T | C | G |
    M          | A |   | T | G | G |           | A | T | C | G |
    Pos          1       3   4   5               9  10  11  12
    Diffpos      2       1   1   4               1   1   1   1
    read 1     | A |   | T | G | G |                                                 POS = 1 CIGAR = '1M1D3M'
    read 2             | T | G | G |           | A |                                 POS = 3 CIGAR = '3M3D1M'

Note that we don't talk about SNPs because they don't change the sequence length and so read POS and CIGAR are very
straightforward.


Misc design choices
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
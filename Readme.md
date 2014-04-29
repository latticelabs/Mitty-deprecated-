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


Files and formats
-----------------

## VCF

`vcf2seq.py` currently handles a specific interpretation of the VCF. In the most liberal interpretation of the VCF the
REF sequence corresponds to bases matching the reference starting at POS and those bases are replaced by the bases found
in ALT. Any number of bases in ALT may match the REF in any place (though this is not very useful to us).

Mitty has a more strict interpretation of the VCF. In Mitty's interpretation, at most, the first base
of REF will match with ALT. All other bases must be different. All variants can be coded in this manner:

Say our original sequence is `ATCGGATC`

                 POS   REF     ALT
    Deletion      1   ATCG     A     -> AGATC
    Deletion      2   TCG      .     -> AGATC     (Though this form is interpreted by vcf2seq it is never produced by mutate.py)
    Deletion      1    A       .     -> TCGGATC
    Insertion     0    .      TTT    -> TTTATCGGATC
    Insertion     1    A      AGGG   -> AGGGTCGGATC
    SNP           1    A       G     -> GTCGATC



Dev notes
=========

Plugin system for simulation models
-----------------------------------

Mitty implements models for variant and read simulation as Python modules located in the Plugins directory. The modules
need to expose a few key functions that `mutate.py` and `reads.py` use to determine variant and read characteristics.
Mitty infers the module name from the name given in the parameter .json file. Mitty comes with some stock models for
variant and read simulation which can be used as example code.

### Random number generators, blocked computation etc.
You will note that the stock plugins accept one or more inputs that serve as seeds for internal random number
generators. These seeds need to be specified in the simulation parameter files and ensure reproducibility of
simulations.

Some of the models use several, independently seeded, random number generators for different variables (e.g.
insertion length and insertion position) to avoid unexpected interactions between such simulated variables. The
algorithms are also designed such that the simulation results are independent of block size.


Algorithms
----------
### Variants, reads and CIGARS
One big goal of Mitty is to serve up realistic test data for bioinformatics algorithms, from aligners to variant
callers. Testing whether a variant caller is correctly working on the simulated data is relatively easy: we simply
compare the variant caller's VCF file with the answer book VCF generated by Mitty. It is, however, slightly more involved
to deduce if an aligner is correctly aligning the simulated reads, and to diagnose how the performance of an aligner is
affecting the accuracy of a variant caller. To this end Mitty has a system to compute the correct read position and
CIGAR for each simulated read. This information is stored in the read's qname string so that it is easily accessible
to diagnostic programs.

In order to generate reads based on a given VCF file and a reference sequence we go through a two step process.
We first generate the mutated sequence (`mut_seq`) along with some other information that encodes the difference between
each base in the `mut_seq` and corresponding positions on the `ref_seq`. We then generate reads from the `mut_seq` using
the sidecar information to compute the correct POS values and CIGAR strings for the reads.

The algorithm is best introduced through a series of examples. In the examples the reference sequence is labelled `R` and
the mutated sequence is labeled `M`. The position information for the read is taken from an array `pos` that accompanies
`M` and the CIGAR is written from information carried by another array `dpos` accompanying `M`.

#### Generating `M`, `pos` and `dpos` arrays from VCF entries

Taking our original example sequence `ATCGGATC` we show examples of SNPs, deletions and insertions as well as
combinations of these variants

SNP at start

    POS REF ALT
    1   A   T

    ref  12345678
         ATCGGATC
    mut  TTCGGATC
    pos  12345678     (Position values correspond to those of the original bases that survive in the mutation)
    dpos 11111111     (Simply, the point by point diff of the pos array. We imagine there to be an imaginary base
                       at the end - with pos 9 - which gives us dpos=1 for the last base)


SNP in middle

    POS REF ALT
    4   G   A

    ref  12345678
         ATCGGATC
    mut  ATCAGATC
    pos  12345678
    dpos 11111111


Insertion at the very beginning of the sequence

    POS REF ALT
    0   .   TC

    ref    12345678
           ATCGGATC
    mut  TCATCGGATC
    pos  1112345678    (For each inserted base, the pos value is copied from the next conserved base)
    dpos 0011111111


Insertion after the first base

    POS REF ALT
    1   A   AA

    ref  1 2345678
         A TCGGATC
    mut  AATCGGATC
    pos  122345678
    dpos 101111111


Insertion at end

    POS REF ALT
    8   C   CTT

    ref  12345678
         ATCGGATC
    mut  ATCGGATCTT
    pos  1234567899   (Note that pos of inserted base is pos of the imaginary last base)
    dpos 1111111100    (Note that the last dpos is 0. Recall, we have a last imaginary base with pos 9)

    * Even though the imaginary base position actually appears in the pos string, no read POS ever sees this - a read will
      either have at least one earlier base matching the reference (which will form the read POS) or there will be no
      matching bases, in which case the read will be unmapped.


Deletion of the first base

    POS REF ALT
    1   A   .

    ref  12345678
         ATCGGATC
    mut   TCGGATC
    pos   2345678
    dpos  1111111


Deletion in the middle

    POS REF ALT
    2   TCG  T

    ref  12345678
         ATCGGATC
    mut  AT  GATC
    pos  12  5678
    dpos 13  1111


Adjacent SNPs

    POS REF ALT
    1   A   T
    2   T   C

    ref  12345678
         ATCGGATC
    mut  TCCGGATC
    pos  12345678
    dpos 11111111


Adjacent Insertions

    POS REF ALT
    0   .   T
    1   A   AT

    ref   1 2345678
          A TCGGATC
    mut  TATTCGGATC
    pos  1122345678
    dpos 0101111111


Insertions and deletion adjacent to each other

    POS REF ALT
    0   .   T
    1   A   .

    ref   12345678
          ATCGGATC
    mut  T TCGGATC
    pos  1 2345678
    dpos 1 1111111

This is a kind of gotcha! It's actually a SNP (A -> T @ pos 1) and as we will find out below our POS and CIGAR computations
will be correct but it raises the point that there is more than one way to describe some mutations. In this case a
SNP is the least complex (and therefore correct) way.

`vcf2seq.py` takes in a VCF file and a reference sequence and generates three files

    mut_seq - the complete, literal, mutated sequence
    pos     - a uint16 array the same size as mut_seq, containing coordinate information that can be used to
              generate the POS value for a correctly aligned read
    dpos    - a uint16 array which is the diff of the pos array and that can be used to generate correct CIGAR
              strings for reads.

`vcf2seq.py` iterates over the `ref_seq` copying over each base into `mut_seq` and copying the base's position in
the reference position into `pos`. `dpos` is set to 1 for such runs of regular copying. When we arrive at the
location of a variant (say it is at `k`), we copy over the ALT, then we advance along `ref_seq` as many bases as there
was in REF (this can be zero if REF='.'). Say this brings us to location `m`. We fill `pos` with [`k`, `m`, `m` ...]
based on the length of ALT. We fill `dpos` with [`m-k`, 0, 0, 0 ... ] based on the length of ALT

#### Generating reads with proper `POS` and `CIGAR`s

Consider some of the examples from above:

SNP at start

    POS REF ALT
    1   A   T

    ref  12345678
         ATCGGATC
    mut  TTCGGATC
    pos  12345678
    dpos 11111111     example reads
          TCGG      POS = 2 CIGAR = '4M'
             GATC   POS = 5 CIGAR = '4M'



Insertion at the very beginning of the sequence

    POS REF ALT
    0   .   TC

    ref    12345678
           ATCGGATC
    mut  TCATCGGATC
    pos  1112345678
    dpos 0011111111   example reads
         TCAT        POS = 1 CIGAR = '2I2M' -> POS is the position of the first matching base
                                               2M comes from counting up the consecutive 1s
                                               2I comes from counting up the consecutive 0s
          CATC       POS = 1 CIGAR = '1I3M'


Adjacent Insertions

    POS REF ALT
    0   .   T
    1   A   AT

    ref   1 2345678
          A TCGGATC
    mut  TATTCGGATC
    pos  1122345678
    dpos 0101111111   example reads
         TATT         POS = 1 CIGAR = '1I1M1I1M' pos is the value corresponding to the first


Insertions equal to or longer than read

    POS REF ALT
    3   C   CATAT

    ref  123    45678
         ATC    GGATC
    mut  ATCATATGGATC
    pos  123444445678
    dpos 111000011111 example reads
            ATAT      Unmapped -> if all dpos values are 0, this is an unmapped read
         ATCA         POS = 1 CIGAR = '3M1I'
               TGGA   POS = 4 CIGAR = '1I3M'


Deletion in the middle

    POS REF ALT
    2   TCG  T

    ref  12345678
         ATCGGATC
    mut  AT  GATC
    pos  12  5678
    dpos 13  1111   example reads
         AT  GA     POS = 1 CIGAR = '2M2D2M' - the 2D comes from 3 - 1. The 1 goes into the Ms (so we get 2M before that)
          T  GAT    POS = 2 CIGAR = '1M2D3M'



Insertions and deletion adjacent to each other (the gotcha)

    POS REF ALT
    0   .   T
    1   A   .

    ref   12345678
          ATCGGATC
    mut  T TCGGATC
    pos  1 2345678
    dpos 1 1111111  example reads
         T TCG      POS = 1 (The first position with dpos > 0) CIGAR = '4M'
                    Note how the CIGAR is correct.


`reads.py` generates simulated reads from `mut_seq` based on the read model. Using the `pos` and `dpos` arrays it
also generates appropriate alignment information that is stored in the qname string as a POS value and a CIGAR string.
(Note that while the BAM specs do not place a limit on the length of the qname string both Tablet and IGV expect a
string with length < 255 characters. It is possible that the qname will exceed this and you won't be able to open a
set of simulated reads using tools that arbitrarily limit the qname). If no `pos` and `dpos` files are supplied
`reads.py` assumes we are taking reads from a reference sequence and the POS values are actual positions of the reads
and all the cigars are of the form `100M` (For e.g. 100 base reads).

Computing POS: Step through the `dpos` array corresponding to the read and set POS corresponding to the first non-zero
value. If all values of `dpos` are zero, it is an unmapped read and POS and CIGAR strings are undefined

Computing the CIGAR: Only do this for mapped reads.
1. Setup three counters, I,D,M and initialize to zero.
2. Step through the relevant part of the `dpos` array.
3. If the value is 1 increment the M counter, flush other counters
4. If the value is 0 increment the I counter and flush other counters
5. If the value is > 0 then increment the M counter, flush it and then increment the D counter by 1 less than the value
6. Continue (2) until we are at the end of the read.
7. At the end, flush the counter that remains.

Flushing a counter involves writing out xM or xD or xI into the cigar string, where `x` is the (non-zero) counter value.


Misc design choices
--------------
### Choice to output a mutated sequence as a whole
Though generating a whole mutated sequence uses a lot of disk space, I chose this approach as it ended up being simpler
than coming up with an algorithm for generating reads on the fly based ona  VCF file. In the future the code may be
converted to do on the fly generation.

### Parameter files for mutate
1. I chose to use parameter files because we often want to rerun experiments and it became clear early on that there would
be a lot of parameters.
1. I chose to use python for the parameter file for parsimony and flexibility
1. The parameter distribution between file and command line was based on predictions of which parameters we could
experiment with most during testing
1. At this time I do not know whether having everything on the commandline would be better for PIPITOR or if param files
are preferred for Platform integration, but either way is a short code reorganization that can be done quickly at the
time of integration.

### POS and DIFFPOS files
These are simple binary files carrying unsigned 4 byte int information. This is enough to handle index/index diff sizes
for human genome sizes, though if we ever work on heavily mutated specimens of the loblolly pine, perhaps we have to
go to 8 byte ints.



Disorganized - don't read below this
====================================



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
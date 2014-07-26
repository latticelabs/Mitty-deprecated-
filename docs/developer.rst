Developer/Advanced user notes
=============================

HDF5 for Whole genome file
--------------------------
Mitty stores sequences from a genome in a whole genome file. This is a standard HDF5_ file with datasets corresponding
to chromosome sequences with associated metadata.

.. _HDF5: http://en.wikipedia.org/wiki/Hierarchical_Data_Format

* The HDF5 file is gzip compressed, so it is the same size as a comparable multi-fasta file,
* Any sequence in the file can be accessed directly, rather than via sequential access [#block]_
* It can hold metadata about the genome and the individual chromosomes in a natural, well defined manner
* Sequences are organized into chromosomes and copies of chromosomes, allowing us to store simulated genomes with special characteristics such as polysomies_.
* Extra data can be stored that is used further down the simulation chain (e.g. POS arrays for generating proper CIGAR and POS values for simulated reads)

.. _polysomies: http://en.wikipedia.org/wiki/Polysomy
.. [#block] Early versions of Mitty used block processing for everything. I changed this to whole chromosome based processing because for Human use cases individual chromosomes are small enough to fit even on modest machines. This considerably simplifies the code and concepts.

Whole genome file structure
---------------------------
The structure of the whole genome file is as follows::

    /sequence/<chrom>.attrs['seq_id'] -> sequence id e.g. NIH accession number
    /sequence/<chrom>/<cpy>  -> sequence data for this chromosome and copy
    /sequence/<chrom>/<cpy>.attrs['reference']  -> boolean. True if this is copied from the reference genome

    pos/<chrom>/<cpy>  -> pos array of u4 (see reads.py for its utility)


mutate.py
---------
Internally, Mitty stores each variant as a group of operations

(type, het, footprint)
()



The reference is assumed to be haploid (we always use copy 1)

Mitty stores the actual generated variant in a systematic format in a hdf5 file. Conceptually, each variant is stored
as follows:

   _ type name
  /
(type, het, footprint, )
        \            \
         \            \_ list of variant description tuples
          \
              0 -> no variant (due to arbitration, see below)
              1 -> variant on copy 1
              2 -> variant on copy 2
              3 -> variant on both copies

Most commonly, there is only one element in the list of variant descriptions and this corresponds exactly with a VCF
file entry. Complicated variants are expressed as a sequence of insertions and deletions.

Each variant model returns a variant file of proposed variants in this format. Mitty then arbitrates between all the
variants to make sure they don't clash.

Variant types:

SNP             - footprint is 1
Insert          - footprint is 0
Delete          - footprint is N
Inversion       - footprint is N
Repeat          - footprint is N
Translocation   - foorprint is N

Mitty writes out variants in VCF format (using the variant2vcf tool).

Each model should have a .variant method that is passed the reference genome data and any general and specific
parameters it needs. Each model's .variant method is called in turn and it fills out the variant data structure
(on an hdf5 file, due to the potential size of the data) and returns it to mutate.py

mutate.py is responsible for arbitrating between variants that clash by using a genome-wide mask.



Passing whole genome file to mutation plugins
---------------------------------------------
This allows plugins to do non-local things, like translocations, which they would not be able to do if they only recieved
a sequence.



Choice to output a mutated sequence as a whole
----------------------------------------------
VCF files (especially with the literal phase information as used by Mitty) form a complete description of a genome (in the form of diffs to the reference, haploid, sequences). However, for the purposes of taking reads the whole mutated genome is written out explicitly. I chose this approach as it ended up being simpler than implementing an algorithm for generating reads on the fly based on a VCF file. In the future the code may be modified to do on the fly generation.


Random number seeds
-------------------
All stock plugins employ explicit random number seeds. Random number generators for different parameters of the simulation are decoupled (independent) from each other and each takes its own seed (Though all plugins can take a single master seed which they use to generate the required number of individual seeds). This is an important design choice that allows exact reproducibility of simulations and allows us to avoid couplings between parameters which might otherwise crop up if the same random number generator was used for all the simulation variables.


vcf2seq - why we store VCF details
----------------------------------
In some use cases [#localreads]_ we need to know which parts of the mutated genome are mutated. For this reason we store the positions of mutations (in `/variants/pos/`) and the details of the variant (`/variants/codes/`).

.. [#localreads] The `reads.py` program's `localreads` option uses this, for example, to generate reads from only the insertions.


Algorithms for computing correct CIGAR and POS values for simulated reads
--------------------------
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
the mutated sequence is labeled `M`. The information for setting the POS and CIGAR for the read is taken from an
array `pos` that accompanies `M`

Generating `pos`
~~~~~~~~~~~~~~~

Let ``R = ACTGACTG``

Consider a single base insertion at position 1::

    POS REF ALT
    1   A   AT

         1 2345678
    R    A CTGACTG
    M    ATCTGACTG
    pos  1223456789


Consider a multiple base insertion at position 1::

    POS REF ALT
    1   A   ATT

         1  2345678
    R    A  CTGACTG
    M    ATTCTGACTG
    pos  12223456789


Consider a multiple base insertion at last position::

    POS REF ALT
    8   G   GTT

         12345678
    R    ACTGACTG
    M    ACTGACTGTT
    pos  12345678999

Consider a multiple base deletion::

    POS REF ALT
    2   CTG  C

         12345678
    R    ACTGACTG
    M    AC  ACTG
    pos  12  56789

Consider a SNP, an insertion and a deletion::

    POS REF ALT
    2   C   T
    4   G   GTT
    6   CTG C

         1234  5678
    R    ACTG  ACTG
    M    ATTGTTAC
    pos  123455569


`pos` is generated by copying over the index from `R`. When we encounter an insertion we copy over the index of the next
reference base as many times as there is an insertion. Deletions are simply skipped. For the purposes of computing `pos`
we also add an imaginary base position at the end of the reference sequence (9 in this case)

Generating CIGARS and POS for reads from `pos`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider our last example and some reads from `M`::

         1234  5678
    R    ACTG  ACTG
    M    ATTGTTAC
    pos  123455569
         ++++---------> POS = 1 (The first pos value we encounter)
                        CIGAR = 4M  (2-1=1 -> 1M
                                     3-2=1 -> 1M
                                     4-3=1 -> 1M
                                     5-4=1 -> 1M)

    M    ATTGTTAC
    pos  123455569
          ++++--------> POS = 2 (The first pos value we encounter)
                        CIGAR = 3M1I  (3-2=1 -> 1M
                                       4-3=1 -> 1M
                                       5-4=1 -> 1M
                                       5-5=0 -> 1I)

    M    ATTGTTAC
    pos  123455569
           ++++-------> POS = 3
                        CIGAR = 2M2I  (4-3=1 -> 1M
                                       5-4=1 -> 1M
                                       5-5=0 -> 1I
                                       5-5=0 -> 1I)

    (A basic aligner would actually soft-clip these last two reads giving us 3M1S and 2M2S)

    M    ATTGTTAC
    pos  123455569
             ++++-----> POS = 5
                        CIGAR = 2I2M  (5-5=0 -> 1I
                                       5-5=0 -> 1I
                                       6-5=1 -> 1M
                                       9-6=3 -> 1M + 2D) The D only comes into play if our read crosses the deletion

To see how a deletion affects our POS and CIGAR consider another previous example::

    POS REF ALT
    2   CTG  C

         12345678
    R    ACTGACTG
    M    AC  ACTG
    pos  12  56789
         ++  ++-------> POS = 1
                        CIGAR = 2M2D2M  (2-1=1 -> 1M
                                         5-2=3 -> 1M + 2D The 2D comes into play because the read crosses the boundary
                                         6-5=1 -> 1M
                                         7-6=1 -> 1M)

Example of an unmapped read::

    POS REF ALT
    2   C  CAATTGG

         12      345678
    R    AC      TGACTG
    M    ACAATTGGTGACTG
    pos  123333333456789
           ++++-------> POS = 3
                        CIGAR = 4I  (3-3=0 -> 1I
                                     3-3=0 -> 1I
                                     3-3=0 -> 1I
                                     3-3=0 -> 1I)
    For a read to be mapped, there has to be at least one M. Since there are no Ms we discard the POS and CIGAR as this
    is an unmapped read

``reads.py`` generates simulated reads from ``mut_seq`` based on the read model. Using the `pos` arrays it
also generates appropriate alignment information (POS and CIGAR) that is stored in the qname string.
(Note that while the BAM specs do not place a limit on the length of the qname string both Tablet and IGV expect a
string with length < 255 characters. It is possible that the qname will exceed this and you won't be able to open a
set of simulated reads using tools that arbitrarily limit the qname). If no `pos` file is supplied `reads.py` assumes
we are taking reads from a reference sequence and the POS values are actual positions of the reads and all the cigars
are of the form `100M` (For e.g. 100 base reads).

Computing POS: For every read, the POS value is simply the index from `pos` corresponding to the first base of the read
EXCEPT for unmapped reads.

Computing the CIGAR:

1. Initialize the base counter to `None`, set mapped flag to `False`
2. Step through the each base of the read and look at the difference in `pos` values `dp`
3. If `dp==1`, if the counter is any thing other than `M`, flush it. Set or increment counter as `M`. Set mapped flag to `True`
4. If `dp==0`, if the counter is other than `I`, flush it. Set or increment counter as `I`
5. If `dp>1`, if the counter is other than `M`, flush it. Set and flush counter as `M`, set counter as `D` to be dp-1
6. Continue from 2 until done.
7. Flush any counter other than `D`
8. If the mapped flag is `False` reset POS and CIGAR - this is an unmapped read.

You can "read along" to these examples by running `python reads.py test -v` and seeing how different functions in
`reads.py` implement these algorithms



### POS files
These are simple binary files carrying unsigned 4 byte int information. This is enough to handle index/index diff sizes
for human genome sizes, though if we ever work on heavily mutated specimens of the loblolly pine, perhaps we have to
go to 8 byte ints ...




Python's native mmap can't do proper offsets ... should we use numpy?

Running tests
-------------

Running all tests::

    nosetests tests  -v


Including the few doctests there are::

    nosetests mitty --with-doctest -v

Including specific doctests::

    nosetests mitty/Plugins/Mutation --with-doctest -v

Running specific tests::

    nosetests tests.vcf2seq_test:test_assemble_sequences_hetero_same_locus_del -v


Generating documentation
------------------------

`sphinx-apidoc mitty/ -o docs` from the root directory

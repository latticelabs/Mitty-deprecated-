Developer/Advanced user notes
=============================

.. toctree::
    :maxdepth: 2

    variants
    genomes



De Novo variant collision detection
-----------------------------------

When Mitty places de novo variants it prevents collisions with earlier variants by using a mask. I ended up going with
bitarray because it gave the best space/time tradeoff.


Speed optimizing denovo.py
--------------------------

Initially I was using scipy.sparse.lil_matrix for the collision mask, because I figured that this was the best way to economize on
space. I found, however, by simple profiling, that lil was taking what I considered to be a lot of time. I tried various other
sparse matrix types but lil was the fastest. I was especially bothered by the overhead of the lil __getitem__ call which
I thought was taking way too much cumulative time.

I was afraid to use regular arrays (regular numpy arrays gave a 5x speed up compared to sparse arrays) because I thought
this would not scale spacewise, because a whole genome simulation would take about 7 GB for the mask alone. Then I tried
bitarray. This turned out to be really fast though it does take some time during initialization.

Going from scipy.sparse to bitarray gave a 7x speed boost, knocking down execution time from 70s to 10s. Very interestingly,
when I replaced the `Variant` data type (which was a `named_tuple`) with a ctypes struct class the total time was knocked
down to 3s. I was shocked that a data type could be so inefficient, but it is used everywhere, and all the time, so
it is logical that any inefficiency here creates such a great effect.

So overall, for sprinkling about 130000 SNPs on chromosome 11, denovo.py went from 70s to 3s by some simple alterations.

Profile results: chr11 130000 SNPs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The following profile runs demonstrate why I ended up using `bitarray` for the collision detection mask

typical command line::

    python -m cProfile -o bitarray_stats.bin ../mitty/denovo.py  --wg ../examples/Human/Out/chr_11.h5  --vcf test.vcf.gz  --param_file ../examples/denovo/snp.json  --master_seed=1 -v
    DEBUG:__main__:Reference file ../examples/Human/Out/chr_11.h5
    DEBUG:mitty.Plugins.variants.snp_plugin:Used master seed to generate seeds 90950894, 39770206, 70599635, 88188254
    DEBUG:__main__:131015 of 131131 variants placed
    DEBUG:variation:Sorting /var/folders/s_/bhgzsrbj7cj220_s1cnpzgvm0000gn/T/tmpyenywa.vcf to test.vcf
    DEBUG:variation:Compressing and indexing test.vcf to test.vcf.gz

sparse matrix::

    In [53]: p = pstats.Stats('sparsearray_stats.bin')
    In [54]: p.strip_dirs().sort_stats('tottime').print_stats(30)
    Thu Aug  7 08:11:20 2014    sparsearray_stats.bin

             49637097 function calls (49636390 primitive calls) in 73.646 seconds

       Ordered by: internal time
       List reduced from 1047 to 30 due to restriction <30>

       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
       655126    7.340    0.000   18.639    0.000 stride_tricks.py:36(broadcast_arrays)
            1    6.791    6.791   65.445   65.445 denovo.py:120(arbitrate_variant_collisions)
       917388    5.115    0.000    8.444    0.000 stride_tricks.py:22(as_strided)
      3407329    4.716    0.000    4.716    0.000 {numpy.core.multiarray.array}
       458694    4.112    0.000   22.274    0.000 sputils.py:329(_index_to_arrays)
       262262    3.974    0.000    3.974    0.000 {scipy.sparse._csparsetools.lil_fancy_get}
       196432    2.690    0.000    2.690    0.000 {scipy.sparse._csparsetools.lil_fancy_set}
       917388    2.436    0.000    3.810    0.000 sputils.py:310(_check_boolean)
     10683970    2.344    0.000    2.644    0.000 {isinstance}
       458694    2.230    0.000    5.354    0.000 sputils.py:244(_unpack_index)
       262262    2.170    0.000   36.775    0.000 lil.py:228(__getitem__)
       196432    1.922    0.000   21.189    0.000 lil.py:270(__setitem__)
            3    1.759    0.586    1.759    0.586 {method 'rand' of 'mtrand.RandomState' objects}
       262263    1.759    0.000    8.985    0.000 lil.py:86(__init__)
      2948597    1.525    0.000    5.519    0.000 numeric.py:392(asarray)
       524526    1.453    0.000    2.280    0.000 fromnumeric.py:2412(ndim)
       524548    1.449    0.000    1.449    0.000 {numpy.core.multiarray.empty}
       458694    1.146    0.000    3.341    0.000 shape_base.py:8(atleast_1d)
       458694    0.956    0.000    0.956    0.000 {scipy.sparse._csparsetools.prepare_index_for_memoryview}
       262263    0.913    0.000    4.245    0.000 sputils.py:204(isshape)
       458694    0.867    0.000    0.867    0.000 {method 'reshape' of 'numpy.ndarray' objects}
       458694    0.815    0.000    0.861    0.000 sputils.py:272(_check_ellipsis)
      2489903    0.806    0.000    1.426    0.000 base.py:861(isspmatrix)
       458694    0.795    0.000    0.795    0.000 {numpy.core.multiarray.arange}
            1    0.789    0.789    0.789    0.789 {posix.waitpid}
    7015738/7015606    0.696    0.000    0.696    0.000 {len}
       655075    0.599    0.000    0.599    0.000 collections.py:54(__setitem__)
       262263    0.575    0.000    0.607    0.000 lil.py:132(set_shape)
            2    0.562    0.281    3.391    1.695 snp_plugin.py:51(variant_generator)
       262262    0.551    0.000    0.690    0.000 lil.py:181(getnnz)

numpy array::

    In [51]: p = pstats.Stats('numpyarray_stats.bin')
    In [52]: p.strip_dirs().sort_stats('tottime').print_stats(30)
    Thu Aug  7 08:08:02 2014    numpyarray_stats.bin

             4222136 function calls (4221429 primitive calls) in 10.668 seconds

       Ordered by: internal time
       List reduced from 1015 to 30 due to restriction <30>

       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            1    2.501    2.501    2.670    2.670 denovo.py:120(arbitrate_variant_collisions)
            3    1.764    0.588    1.764    0.588 {method 'rand' of 'mtrand.RandomState' objects}
            1    0.787    0.787    0.787    0.787 {posix.waitpid}
            2    0.575    0.288    3.407    1.704 snp_plugin.py:51(variant_generator)
       655075    0.566    0.000    0.566    0.000 collections.py:54(__setitem__)
            3    0.536    0.179    0.536    0.179 {method 'nonzero' of 'numpy.ndarray' objects}
       131015    0.509    0.000    1.550    0.000 _abcoll.py:526(update)
       131015    0.406    0.000    2.099    0.000 collections.py:38(__init__)
            1    0.364    0.364    3.377    3.377 variation.py:77(vcf_save)
            1    0.335    0.335    0.335    0.335 {method 'read' of 'h5py.h5d.DatasetID' objects}
       131015    0.221    0.000    2.375    0.000 <string>:24(_asdict)
       131193    0.202    0.000    0.202    0.000 {method 'format' of 'str' objects}
       131015    0.180    0.000    0.287    0.000 abc.py:128(__instancecheck__)
       262262    0.169    0.000    0.169    0.000 {numpy.core.multiarray.count_nonzero}
       393058    0.143    0.000    0.143    0.000 _weakrefset.py:70(__contains__)
            1    0.139    0.139    0.139    0.139 {pysam.ctabix.tabix_compress}
       131029    0.128    0.000    0.128    0.000 {map}
       262157    0.126    0.000    0.126    0.000 {built-in method __new__ of type object at 0x100183140}
       131020    0.121    0.000    0.121    0.000 {zip}
       131015    0.111    0.000    0.357    0.000 <string>:28(_replace)
       131057    0.084    0.000    0.084    0.000 {hasattr}
       131021    0.079    0.000    0.079    0.000 {method 'write' of 'file' objects}
    131027/131026    0.076    0.000    0.133    0.000 abc.py:148(__subclasscheck__)
       132944    0.071    0.000    0.358    0.000 {isinstance}
       131015    0.064    0.000    0.117    0.000 <string>:12(_make)
          2/1    0.061    0.030   10.567   10.567 denovo.py:69(<module>)
            1    0.050    0.050    0.050    0.050 {pysam.ctabix.tabix_index}
    528191/528059    0.042    0.000    0.042    0.000 {len}
       131131    0.040    0.000    0.125    0.000 <string>:8(__new__)
       131530    0.020    0.000    0.020    0.000 {getattr}



bitarray::

    In [55]: p = pstats.Stats('bitarray_stats.bin')
    In [56]: p.strip_dirs().sort_stats('tottime').print_stats(30)
    Thu Aug  7 08:16:50 2014    bitarray_stats.bin

             4222178 function calls (4221471 primitive calls) in 10.629 seconds

       Ordered by: internal time
       List reduced from 1015 to 30 due to restriction <30>

       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            3    1.785    0.595    1.785    0.595 {method 'rand' of 'mtrand.RandomState' objects}
            1    1.775    1.775    1.793    1.793 denovo.py:123(arbitrate_variant_collisions)
            1    0.787    0.787    0.787    0.787 {posix.waitpid}
            1    0.662    0.662    0.662    0.662 denovo.py:85(initialize_mask)
       655075    0.593    0.000    0.593    0.000 collections.py:54(__setitem__)
            2    0.578    0.289    3.406    1.703 snp_plugin.py:51(variant_generator)
       131015    0.533    0.000    1.622    0.000 _abcoll.py:526(update)
            3    0.524    0.175    0.524    0.175 {method 'nonzero' of 'numpy.ndarray' objects}
       131015    0.438    0.000    2.206    0.000 collections.py:38(__init__)
            1    0.387    0.387    3.561    3.561 variation.py:77(vcf_save)
            1    0.336    0.336    0.336    0.336 {method 'read' of 'h5py.h5d.DatasetID' objects}
       131015    0.236    0.000    2.500    0.000 <string>:24(_asdict)
       131194    0.218    0.000    0.218    0.000 {method 'format' of 'str' objects}
       131015    0.188    0.000    0.298    0.000 abc.py:128(__instancecheck__)
       393058    0.147    0.000    0.147    0.000 _weakrefset.py:70(__contains__)
            1    0.145    0.145    0.145    0.145 {pysam.ctabix.tabix_compress}
       131029    0.135    0.000    0.135    0.000 {map}
       131015    0.119    0.000    0.374    0.000 <string>:28(_replace)
       131020    0.119    0.000    0.119    0.000 {zip}
       262157    0.118    0.000    0.118    0.000 {built-in method __new__ of type object at 0x100183140}
       131057    0.089    0.000    0.089    0.000 {hasattr}
       131021    0.082    0.000    0.082    0.000 {method 'write' of 'file' objects}
    131027/131026    0.078    0.000    0.137    0.000 abc.py:148(__subclasscheck__)
       132950    0.075    0.000    0.373    0.000 {isinstance}
       131015    0.067    0.000    0.120    0.000 <string>:12(_make)
          2/1    0.062    0.031   10.532   10.532 denovo.py:69(<module>)
            1    0.052    0.052    0.052    0.052 {pysam.ctabix.tabix_index}
    528191/528059    0.042    0.000    0.042    0.000 {len}
       131131    0.039    0.000    0.116    0.000 <string>:8(__new__)
       131530    0.022    0.000    0.022    0.000 {getattr}


bitarray with ctypes class for Variation::

    In [45]: p = pstats.Stats('bitarray_stats.bin')
    In [46]: p.strip_dirs().sort_stats('tottime').print_stats(30)
    Fri Aug  8 06:18:49 2014    bitarray_stats.bin

             420785 function calls (420089 primitive calls) in 3.633 seconds

       Ordered by: internal time
       List reduced from 923 to 30 due to restriction <30>

       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            1    0.884    0.884    0.884    0.884 {posix.waitpid}
            1    0.657    0.657    0.659    0.659 denovo.py:84(initialize_mask)
            2    0.594    0.297    1.098    0.549 snp_plugin.py:51(variant_generator)
            1    0.415    0.415    0.415    0.415 {method 'read' of 'h5py.h5d.DatasetID' objects}
            1    0.289    0.289    0.306    0.306 denovo.py:114(arbitrate_variant_collisions)
            1    0.178    0.178    0.235    0.235 variation.py:103(vcf_save)
            1    0.147    0.147    0.147    0.147 {pysam.ctabix.tabix_compress}
            1    0.069    0.069    0.069    0.069 {pysam.ctabix.tabix_index}
            5    0.064    0.013    0.064    0.013 {zip}
       131115    0.057    0.000    0.057    0.000 {method 'write' of 'file' objects}
          2/1    0.055    0.027    3.541    3.541 denovo.py:69(<module>)
            3    0.021    0.007    0.021    0.007 {open}
       262218    0.017    0.000    0.017    0.000 {method 'any' of 'bitarray._bitarray' objects}
            1    0.015    0.015    0.015    0.015 {method 'poisson' of 'mtrand.RandomState' objects}
            3    0.014    0.005    0.087    0.029 __init__.py:10(<module>)
            1    0.010    0.010    1.415    1.415 denovo.py:150(add_variant_model_to_genome)
            9    0.010    0.001    0.058    0.006 __init__.py:1(<module>)
            2    0.006    0.003    0.006    0.003 {method 'read' of 'file' objects}
            1    0.005    0.005    0.006    0.006 __init__.py:88(<module>)
            2    0.005    0.002    0.006    0.003 {__import__}
            1    0.004    0.004    3.451    3.451 denovo.py:195(main)
          273    0.003    0.000    0.004    0.000 function_base.py:2945(add_newdoc)
            2    0.003    0.002    0.034    0.017 variation.py:6(<module>)
            1    0.003    0.003    0.004    0.004 polynomial.py:55(<module>)
            5    0.003    0.001    0.004    0.001 collections.py:288(namedtuple)
            7    0.003    0.000    0.003    0.000 {method 'sub' of '_sre.SRE_Pattern' objects}
       114/42    0.003    0.000    0.007    0.000 sre_parse.py:379(_parse)
            1    0.002    0.002    0.003    0.003 hermite.py:59(<module>)
            1    0.002    0.002    0.007    0.007 util.py:11(het)
            1    0.002    0.002    0.003    0.003 laguerre.py:59(<module>)



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


JIT expansion of variant sequences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Some read models have a fixed


### POS files
These are simple binary files carrying unsigned 4 byte int information. This is enough to handle index/index diff sizes
for human genome sizes, though if we ever work on heavily mutated specimens of the loblolly pine, perhaps we have to
go to 8 byte ints ...




Python's native mmap can't do proper offsets ... should we use numpy?

Running tests
-------------

Running all tests::

    nosetests -v


Including the few doctests there are::

    nosetests mitty --with-doctest -v

Including specific doctests::

    nosetests mitty/Plugins/Mutation --with-doctest -v

Running specific tests::

    nosetests tests.vcf2seq_test:test_assemble_sequences_hetero_same_locus_del -v

Running with coverage::

    nosetests --with-coverage --cover-package=mitty


Profiling the code
------------------

::

    python -m cProfile -o denovostats.bin ../mitty/denovo.py  --wg ../examples/mutate/Out/chimera.h5  --vcf test.vcf.gz  --param_file ../examples/denovo/snp.json -v --master_seed=1

    import pstats
    p = pstats.Stats('denovostats.bin')
    p.strip_dirs().sort_stats('cumulative').print_stats(10)


Generating documentation
------------------------

`sphinx-apidoc mitty/ -o docs` from the root directory

Roadmap
-------

v 4.0.0 April 2016
  - Fuzzer integration


v 3.0.0 Jan 2016
  - Preliminary integration with platform

    - Minimum: should be able to run on platform and produce reports
    - extensions to CWL?
    - should be able to interpret pipeline description to allow us to build


v 2.0.0 Dec 2015
  - Variant Call comparator

    - Includes report generator
    - modular: should support several types of comparators
    - visualizations for comparisons
    - meta analysis of quantified results

  - Haplotype block generator
  - Better simulation management on platform via metadata(?)

Planned changes
---------------
* in perfectbam, use the mismatch count field to produce a summary of #mismatches in a read vs number of reads
* setup.py: figure out how to make "extras_require" work.
* Alindel should be able to split indel accuracy by known vs novel and save those separately
* Alindel should be able to split indel accuracy based on reads under graph variants but not bearing sample variants.
* All plugins should provide commandline stubs that can generate parameter file fragments based on command line
  inputs. This is for compatibility with CWL (?)
* Brandi's use case: compare pipelines against each other.
* Benchmarking: augment circle and matrix plot to plot sub-regions, with the possibility to add
  positions of variants. This is needed for ad hoc questions related to whether there are sample variants
  or graph variants at certain positions
* Reads parameter file should take explicit file names (?)
* change CLIs to allow parameter files to be fed in from std in (piped as a file)
* vcf-write should have option to remove unzipped .vcf file
* variants only reads - need to properly handle end of variants
* finish implementing ad hoc filtering for variants (standard population plugin)

   - het/hom filtering in standard population model
   - range filter for variations (put as core spec - like for reads - rather than in population model?)

Correctness
-----------
* Need to generate reads more randomly. Can perhaps parallelize generation?

   - For every block, generate reads from all chroms and then, when writing them out, pick each read randomly from
     a chrom
   - Run several such processes in parallel, generating separate fastq files. In the end, concatenate them together
     Works for both gzipped and unzipped the same!



Efficiency changes
------------------
* Load chroms one at a time
* Drop reference base(s) from variant list, store SNPs separately (will HDF5 optimize for this?)
* Parallelize all operations by Chrom - will HDF5 allow parallel writes to separate data sets?. May be switch to having
  spearate files for each chrom.




Possible changes
----------------
* Interpret lower case in FASTA as repeat regions and have option to avoid these for variants and reads
* (non-backward compatible) change reads.json format such that read model name is key for model parameters
  much like how it is done for variant models.
* Heuristic for coverage per block?

Partial Changelog
-----------------
(Please see git commit log for detailed commentary)

**1.39.0.dev0**
  * `genome-file` summary command now can give variant counts of multiple samples in a table

**1.38.2.dev0**
  * New sub command for genomes to convert VCF into genome DB

**1.35.0.dev0**
  * Using an efficient filter to discard deletions that contain 'N's anywhere

**1.34.0.dev0**
  * Overhauled genome DB data set organization
    (HDF5 file org is now different and breaks compatibility with earlier versions)

**1.30.0.dev0**
  * genomes and reads modifed so that I/O files can be overridden from the command line.


**1.29.0.dev0**

2015.11.11
  * Cythonized bottlenecks in genome generation

**1.27.0.dev0**

2015.11.06
  * Alindel Plot can now infer indel range from data
  * Auto scale lines/circles in misalignment plots
  * alindel_plot should handle case where there are no indels (log scaling fails)


**1.26.1.dev0**

2015.11.05
  * Bugfix: creed.read_analyze now properly handles position checking of reads with all I or all S

2015.11.04
  * [Wrappers] Use metadata to keep track of files from different aligner versions
  * [Wrappers] Have perfectbam and alindel and alindel plot operate on lists (doing scatter gather possibly)
  * moved wrapper code into separate project


**1.26.0.dev0**

2015.10.30
  * Variant count from indel analysis now only counts variants with at least one read covering them. This takes care of
the counting problems when we take reads from only one chromosome, or only part of a chromosome etc.
  * Alindel plot now shows pairwise differences in additional panel

**1.25.0.dev0**

2015.10.26
  * Combine multiple (or at least two) BADBAMs to perform intersection and difference analyses. Interactive tool?


2015.10.20
  * matrix plot should show light gray dots for grid points


2015.10.19
  * Implemented option to filter multiple allele loci.


2015.10.14
  * Update plot_align (diff ways to plot mis-alignments) to work with BAM+tags way of saving misalignments

2015.10.11
  * Improved documentation

2015.10.07
  * In read simulator/plugins 'SSS..' for the sequence/phred score strings has been changed into 'O' ('object', like for variants)
  * Some of the read plugin code has been abstracted into a base class, allowing us a standard dtype for the numpy arrays
    and one common helper function (get_zero_reads)

----

**1.18.0.dev0**

2015.10.06
  * Enhancement: Full chain upto indel accuracy plot now works
  * Enhancement: Ad hoc post filters implemented in standard population model.
    het/hom filters still need to be implemented
  * Bugfix: Now have a function return empty read array. This fixes an issue with read array concatenation: If we asked for
    reads from variants only, but there were no variants, we would try to concatenate an empty list which would lead to
    an error. This also fixes the problem that in such a condition the paired-endedness of the file would be uncertain.

----

**1.16.0.dev0**

2015.10.05
  * Read length information added to qname

2015.10.01
  * Instead of making several different files write out the alignment accuracy in the original BAM itself.
    Still produce a perfect BAM as needed

2015.09.29
  * Modified read simulator to allow reads to be generated over a sub-region of a chromosome.
    Coverage is correct. Sub-regions have to be set chromosome-by-chromosome.
    Parameter file format change is backwards compatible. Existing parameter files will work correctly with new version
  * Added flag in read simulator to write gzipped fasta file.
    Existing parameter files will work correctly with new version
Standing tasks
==============
1. Update documentation and test as you read
1. Wrap latest version of mitty for IGOR
    - wrapper should detect .fasta file and converta it if needed
    - wrapper should detect vcf file and compress/index it if needed
    - same for bam files
    - Have a overall wrapper, that allows us to go from sequence to reads in one block


List of Mitty tasks in order of priority
========================================
Definitely email kaushik.ghose@sbgenomics.com if you have a suggestion/feature request to be added

1. Efficiency: make roll cigar cython accept list of reads - will save on the function call overhead!
1. Mutate flag move to chromosome rather than sequence?
1. RNGs - option to autogen seeds from single master seed
1. all programs - move I/O (file names) to command line

DONE
====
1. Create lock for delete length to make sure it does not exceed sequence length (DONE: 7/16/2014)
1. Stock insert plugin - have base frequency match nature? (DONE: 7/1/2014)
1. Test het reads with inserts/deletes and see if the chain makes sense. (DONE: 7/1/2014)
1. Modify reads.py and cheata.py to handle full genome (multiple chromosomes).
    - let reads.py access a list of .smalla files corresponding to chromosomes
    - in the reads qname allow for a sequence id to identify which chromosome the read came from
    - in cheata in the written BAM file make sure to copy over the sequence list
    - as part of the input parameters we should have a json file with seq ids and sequence file names which should
  filter through to the BAM files we produce so that cheata produces a perfectly aligned BAM file (minus long inserts)
  for the whole genome (or as much of it as we care to do) (DONE: 7/1/2014)
1. Heterozygous reads clean up and draw out architecture for whole genome stuff. Start with Readme.md and then
   proceed (DONE: 7/1/2014)
1. Handle the CIGAR business (DONE: 5/2/2014)
1. Have a document catchup marathon (DONE: 5/5/2014)
1. Heterozygous reads (DONE: 5/9/2014)
1. Rewrite mutations so we can physically chain them? This is only when I think about the platform. I like things the
way they are logically and computationwise ... NO! What we do is we have little widgets that generate the required .json file that we then pass into mutate!
Yes! (DONE a long time ago)

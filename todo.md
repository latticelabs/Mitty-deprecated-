List of Mitty tasks in order of priority
----------------------------------------
Definitely email kaushik.ghose@sbgenomics.com if you have a suggestion/feature request to be added

1. Handle the CIGAR business (DONE: 5/2/2014)

1. Have a document catchup marathon

1. Wrap mitty for IGOR

1. Modify reads.py and cheata.py to handle full genome (multiple chromosomes).
    - let reads.py access a list of .smalla files corresponding to chromosomes
    - in the reads qname allow for a sequence id to identify which chromosome the read came from
    - in cheata in the written BAM file make sure to copy over the sequence list
    - as part of the input parameters we should have a json file with seq ids and sequence file names which should
  filter through to the BAM files we produce so that cheata produces a perfectly aligned BAM file (minus long inserts)
  for the whole genome (or as much of it as we care to do)

1. Heterozygous reads
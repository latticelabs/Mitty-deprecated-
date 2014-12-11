Files and conventions
=====================

Reference genome
----------------
Mitty expects the reference genome to be stored as a set of .fasta files under a common directory. Mitty requires each
.fasta file to contain just two lines. The first line is the sequence id and the second line is the sequence itself.
Mitty requires any lower case letters (sometimes used to indicate repeats) to uppercase. See ``examples/data``. A
command line script ``splitta`` is provided that will do this for you.

Variant genome
--------------
Variant genomes are expressed as exact VCF files. These files follow a subset of the VCF specification and are very
strict:

* They only allow SNPs, insertions and deletions
* REF and ALT strings are exact
* The REF must contain at least one reference base.
  Violating this rule will cause off-by-one errors in special cases (see developer docs)
* There is only one sample per file.
* Their phasing information (GT, genotype column) is exact: 1/0 means having a variant on copy 0 of the chromosome and 0/1 means having a variant on copy 1 of the chromosome.
* There are no GT=0/0 entries

Reads
-----
Reads are produced as standard `.fastq` files using sanger `PHRED` scores. The `qname` filed is hijacked, however, to
store the correct alignment answer and looks like this::

    chrom:copy|rN|D|POS1|CIGAR1               unpaired reads

    chrom:copy|rN|D|POS1|CIGAR1|POS2|CIGAR2   paired reads

This is the information that the perfect aligner uses to generate perfect alignments of the reads. The `D` field carries
either a `>` or a `<` and indicates if the read (first read if paired) is forward or reverse complemented respectively.
This is important for understanding how to interpret the sequence string.

The CIGAR string is standard except for the case of reads taken from inside an insertion past the anchor points. Such reads should appear as
*unmapped* when using a linear reference. The CIGAR for such reads take the form::

    >OFFSET

Where OFFSET is the offset of the read from the start of the insertion. The POS value indicates the location of the
insertion (The POS value in the VCF file). In this case if the VCF file uses a non-conventional way of expressing insertions::

    REF   ALT
    .     <insertion>

the location of the insertion (POS value) will be off by one.
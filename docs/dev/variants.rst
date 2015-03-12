.. _variation_struct:

Variations, Chromosomes and Genomes
===================================
Internally a variant is represented by a simple structure called :py:class:`mitty.lib.variation.Variant`.
It is implemented as a Cython class as we need to create many instances of this structure. All the variants
relevant to a given simulation are stored in a **master list**. The master list is a Python dictionary indexed by the
position of the variant + which number allele it is. For details please see the code. Most users should not need
to know the details of the indexing.

Biologically, a given variant may be shared amongst many individual samples. What differs from sample to sample is
the copy of the chromosome the variation is found on (the genotype). Mirroring this, we have
:py:class:`mitty.lib.variation.VariantInstance` which carries a pointer to a variation in the master list and genotype
data. A chromosome from a sample is a linked list of ``VariantInstance`` objects maintained in position order called a
``Sample`` (:py:class:`mitty.lib.variation.Sample`).

This split in data, unimportant when we are dealing with a single sample in isolation, can save us a lot of memory and
computation speed when we are dealing with large population simulations. Without this separation in implementation,
each time we create a sample we would waste time and space creating identical copies of variants.
Using the current implementation, a variant in a sample is a four byte index and one byte of genotype data
Given the bottle necked nature of human variation, this is an advantageous architecture for the program.

A genome is simply a dictionary of ``Sample``s with a corresponding dictionary of master lists the individual
``SampleVariants`` point to.
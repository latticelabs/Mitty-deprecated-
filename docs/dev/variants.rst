.. _variation_struct:

Variations
==========
Internally variants are represented by a simple structure called :py:class:`mitty.lib.variation.Variation`.
The code is implemented as a cython class as we need to create many instances of this structure and we need to
repeatedly copy it from one place to another.

Biologically a given variant may be shared amongst many individual samples. What differs from sample to sample is
the copy of the chromosome the variation is found on and (potentially) its recessive nature and contribution to
phenotype fitness.

Mirroring this, ``Variation`` itself is a wrapper that contains :py:class:`mitty.lib.variation.VariationData` and some
metadata. ``VariationData`` carries information about the location and structure of the variant (POS, REF, ALT) while
the metadata stores which copy if the chromosome the variation is on, whether it is recessive and its simulated fitness
value.

This distinction, unimportant when we are dealing with a single sample in isolation, can save us a lot of memory and
computation speed when we are dealing with large population simulations. Without this separation in implementation,
each time we create a sample we would waste time and space creating identical copies of variants.

Using the current implementation, a variant in a sample is three bytes of metadata indicating which copy(ies) of the
chromosome the variant is on, whether it is recessive and a fitness value. There is then a pointer to the actual (heavy)
variant data with (POS, REF and ALT).

Given the bottle necked nature of human variation, this is an advantageous architecture for the program.

Chromosomes and Genomes
-----------------------
A chromosome is represented as a list of Variations. Genomes are stored as a dictionary with chromosome
numbers (1,2,3...) as keys and lists of variants as values.

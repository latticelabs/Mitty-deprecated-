.. _pop_sims:

Population simulations
======================
Variations observed in human populations are due to shuffling of genes during sexual reproduction and (a little) due
to errors introduced during the copying of DNA during meiosis. Mitty can generate a flotilla of VCF files simulating
samples in a population that observe required statistics.

Concepts
--------

Founder populations
~~~~~~~~~~~~~~~~~~~
Population simulations start with a founder population. To create a founder population we construct a dense list of
variants called *common ancestral variants - AV*. We then create 'haplotypes' by


For each sample of the population we take variants, with replacement,
from *AV* and then, on top of that, add freshly generated denovo variants (*dv*) which are unlikely to overlap with *AV*.
A mix of about 95% variants from *AV* and 5% *dv* will result in gross variant statistics that match that reported by
the 1000G project.


Note that this population will not stand up to deeper scrutiny as it will not exhibit any
recombination hot-spots or any meaningful haplotypes.


Generated Populations
~~~~~~~~~~~~~~~~~~~~~
For some requirements, the initial founder population is enough, while for others, we need to simulate sexual
reproduction over several generations. Sexual reproduction mimics mendelian shuffling of variants, recombination
and denovo mutations. A simulation of fitness and recessiveness is also incorporated which can be used to simulates the
weeding out of deleterious mutations. In all stages we keep a log of the transformations performed such that we can
trace the ancestry of an individual genome. These transformations result in genome samples that have haplotypes,
have discoverable hotspots and are amenable to linkage disequilibrium analysis.


Algorithms
----------
(Also see :ref:`denovo_variant_generation`)

The founder



References
----------
.. [myers2006] Myers S, Spencer CC, Auton A, Bottolo L, Freeman C, Donnelly P, McVean G. The distribution and causes of meiotic recombination in the human genome. Biochem Soc Trans. 2006 Aug;34(Pt 4):526-30. doi:10.1042/BST0340526


p1    .     .
p2       .
p3   .  .
p4     .   .

p1 + p2  .  .  .

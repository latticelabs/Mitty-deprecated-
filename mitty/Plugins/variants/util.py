"""Some functions commonly used by mutation plugins are refactored out here for convenience."""
import numpy
from mitty.variation import HOMOZYGOUS, HET1, HET2


def initialize_rngs(*seeds):
  """A pure convenience function that, when passed a set of seeds returns us a bunch of initialized RNGs."""
  return [numpy.random.RandomState(seed=seed) for seed in seeds]


def het(num_vars=0, phet=0.5, het_rng=None, copy_rng=None):
  """This function determines heterozygosity of variants.

  Parameters
  ----------
  num_vars  : int
              How many variants
  phet      : float (0.0 <= phet <= 1.0)
              Probability that a mutation is going to be heterozygous
  het_rng   : object
              Random number generator for determining heterozygosity e.g. numpy.random.RandomState(seed=1)
  copy_rng  : object
              Random number generator for determining which copy the variant is put in
              e.g. numpy.random.RandomState(seed=1)

  Returns
  -------
  het_type  : int array
              #             0      1      2      3
              gt_string = ['0/0', '0/1', '1/0', '1/1']  # The types of genotypes

  Examples
  --------
  >>> het(num_vars=10, het_rng=numpy.random.RandomState(seed=1), copy_rng=numpy.random.RandomState(seed=2))
  array([2, 3, 2, 1, 2, 2, 2, 2, 1, 3], dtype=uint8)
  """
  het_type = numpy.empty((num_vars,), dtype='u1')
  het_type.fill(HOMOZYGOUS)  # Homozygous
  idx_het, = numpy.nonzero(het_rng.rand(het_type.size) < phet)  # Heterozygous locii
  het_type[idx_het] = HET1  # On copy 1
  het_type[idx_het[numpy.nonzero(copy_rng.rand(idx_het.size) < 0.5)[0]]] = HET2  # On copy 2
  return het_type

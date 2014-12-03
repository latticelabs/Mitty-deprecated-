import numpy
from mitty.lib import SEED_MAX

__author__ = 'kghose'


def initialize_rngs(master_seed, n_rngs=4):
  """Return n_rngs initialized from the master_seed"""
  return [numpy.random.RandomState(seed=seed)
          for seed in numpy.random.RandomState(seed=master_seed).randint(SEED_MAX, size=n_rngs)]


def place_poisson(rng, p, end_x):
  """Given a random number generator, a probability and an end point, generate poisson distributed events. For short
  end_p this may, by chance, generate fewer locations that normal"""
  if p == 0.0:
    return numpy.array([])
  est_block_size = end_x * p * 1.2
  these_locs = rng.poisson(lam=1./p, size=est_block_size).cumsum()
  return these_locs[:numpy.searchsorted(these_locs, end_x)]
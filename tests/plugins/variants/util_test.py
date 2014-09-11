import numpy
from numpy.testing import assert_array_equal
from mitty.plugins.variants import util


def place_poisson_test():
  rng = numpy.random.RandomState(seed=1)
  p = .01
  end_p = 10000
  correct_locs = rng.poisson(lam=1./p, size=end_p).cumsum()
  idx, = numpy.nonzero(correct_locs >= end_p)
  correct_locs = correct_locs[:idx[0]]

  rng = numpy.random.RandomState(seed=1)
  computed_locs = util.place_poisson(rng, p, end_p)

  assert_array_equal(correct_locs, computed_locs)
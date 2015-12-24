"""This deletion model generates uniformly distributed deletion lengths for testing."""
import itertools

import numpy as np

import mitty.lib
import mitty.lib.util as mutil
from mitty.plugins.variants import scale_probability_and_validate

import logging
logger = logging.getLogger(__name__)


__example_param_text__ = """
{
  "p": 0.01,           # probability that the deletion will happen at any given base
  "min_len": 10,   # Lower bound on deletion lengths
  "max_len": 1000  # upper bound on deletion lengths
}
"""

_description = __doc__ + __example_param_text__

_example_params = eval(__example_param_text__)


class Model:
  def __init__(self, p=0.01, min_len=1, max_len=1000, ref=None):
    assert 0 <= p <= 1.0, "Probability out of range"
    assert 0 < min_len < max_len, "Check your min_len and max_len definitions"
    self.p, self.del_len_min, self.del_len_max = p, min_len, max_len

  def get_variants(self, ref, p=None, f=None, seed=1, **kwargs):
    """This function is called by the simulator to obtain variants.

    :param ref: reference sequence as a string
    :param chrom: chromosome number (1,2,3,4...)
    :param p: array/list of probability values
    :param f: array/list of frequency values
    :param seed: seed for the random number generators
    :return: 5 arrays/lists/iterables all of the same length
              pos   - position of SNPs
              stop  - stop locations, (pos + 1 for SNPs)
              ref   - reference base,
              alt   - alt base,
              p     - probability value for this variant. These are uniformly distributed random values
    """
    assert 0 < seed < mitty.lib.SEED_MAX
    logger.debug('Master seed: {:d}'.format(seed))

    base_loc_rng, del_len_rng = mutil.initialize_rngs(seed, 2)

    p_eff = scale_probability_and_validate(self.p, p, f)
    del_locs = mutil.place_poisson_seq(base_loc_rng, p_eff, 0, len(ref), ref)
    del_lens = del_len_rng.randint(low=self.del_len_min, high=self.del_len_max, size=del_locs.shape[0])

    del_locs, del_ends, refs, alts = mutil.discard_deletions_in_illegal_regions(ref, del_locs,
                                                                                (del_locs + del_lens + 1).astype('i4'))
    if len(del_locs):
      del_lens = del_ends - del_locs - 1
      p = 0.5 * np.ones(del_lens.size, dtype=float)

    return del_locs, del_ends, refs, alts, p


def test0():
  """Edge case - no variants generated"""
  ref_seq = 'ACTGACTGACTGACTGACTGACTGACTGACTGACTG'
  m = Model(p=0.00001)
  pos, stop, ref, alt, p = m.get_variants(ref_seq, seed=10)
  assert len(pos) == 0  # This should just run and not crash


def test1():
  """Basic test"""
  ref_seq = 'ACTGACTGACTGACTGACTGACTGACTGACTGACTG'
  m = Model(p=0.1)
  pos, stop, ref, alt, p = m.get_variants(ref_seq, seed=10)
  for p, r in zip(pos, alt):
    assert r == ref_seq[p]


if __name__ == "__main__":
  print _description
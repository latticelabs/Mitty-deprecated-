"""This is the stock deletion generator.

Note: This never generates a deletion at the first base of a sequence.

"""
import numpy
import logging
from mitty.lib import SEED_MAX
from mitty.lib.variation import Variation
import mitty.plugins.variants.util as util

logger = logging.getLogger(__name__)


__example_param_text__ = """
{
  "chromosome": [1],    # List of chromosomes to apply the variant to
  "p": 0.01,            # probability that the deletion will happen at any given base
  "phet": 0.5,          # probability that the variant will be heterozygous
  "del_len_lo": 10,      # lower bound on deletion lengths
  "del_len_hi": 20,     # upper bound on deletion lengths
}
"""

_description = """
This is the stock SNP plugin. A typical parameter set resembles
""" + __example_param_text__

_example_params = eval(__example_param_text__)


def variant_generator(ref={},
                      chromosome=None,
                      p=0.01,
                      phet=0.5,
                      del_len_lo=10,
                      del_len_hi=100,
                      master_seed=1,
                      **kwargs):
  assert 0 <= p <= 1.0, "Probability out of range"
  assert 0 <= phet <= 1.0, "Probability out of range"
  assert 0 < del_len_lo <= del_len_hi, "Check deletion length values"
  logger.debug('Master seed: {:d}'.format(master_seed))
  del_loc_rng, del_len_rng, het_rng, copy_rng = util.initialize_rngs(master_seed, 4)

  for chrom in chromosome:
    ref_seq = ref[chrom]  # Very cheap operation
    del_locs = util.place_poisson(del_loc_rng, p, len(ref_seq))
    del_lens = del_len_rng.randint(low=del_len_lo, high=del_len_hi+1, size=del_locs.size)
    het_type = util.het(del_locs.size, phet, het_rng, copy_rng)

    yield {chrom: [Variation(pos + 1, pos + del_len + 1, ref_seq[pos:pos + del_len + 1], ref_seq[pos], het)
                   for het, pos, del_len in zip(het_type, del_locs, del_lens) if ('N' not in ref_seq[pos:pos + del_len + 1])
                   and pos + del_len < len(ref_seq)]}


def test():
  """Should produce variants."""
  ref = {
    1: 'ACTGACTGACTG',
    2: 'ACTGACTGACTGACTGACTGACTG'
  }
  vg = variant_generator(ref, chromosome=[1, 2], p=1.0, del_len_lo=1, del_len_hi=2)
  for v_list in vg:
    assert isinstance(v_list.values()[0][0], Variation), v_list


def test2():
  """The REFs should match the reference seq."""
  ref = {
    1: 'ACTGACTGACTG',
    2: 'ACTGACTGACTGACTGACTGACTG'
  }
  vg = variant_generator(ref, chromosome=[1, 2], p=1.0, del_len_lo=1, del_len_hi=2)
  for v_list in vg:
    assert isinstance(v_list.values()[0][0], Variation), v_list


if __name__ == "__main__":
  print _description
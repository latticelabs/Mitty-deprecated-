"""This is the stock insertion generator.

Note: This never generates a deletion at the first base of a sequence.

"""
__example_param_text = """
{
  "chromosome": [1],
  "phet": 0.5,
  "p": 0.001,
  "ins_len_lo": 5,
  "ins_len_hi": 10
}
"""

_description = """
This insertion model is designed more for testing algorithms than for generating a natural distribution of insertions.
If you need a exponential/power law like distribution of insertions see exp_inserts.
A typical parameter set resembles
""" + __example_param_text

_example_params = eval(__example_param_text)


import numpy
import mitty.plugins.variants.util as util
from mitty.lib.variation import Variation
from mitty.lib import SEED_MAX
import logging
logger = logging.getLogger(__name__)


def variant_generator(ref={},
                      chromosome=None,
                      p=0.01,
                      phet=0.5,
                      ins_len_lo=100,
                      ins_len_hi=10000,
                      master_seed=1,
                      **kwargs):
  assert 0 <= p <= 1.0, "Probability out of range"
  assert 0 <= phet <= 1.0, "Probability out of range"
  logger.debug('Master seed: {:d}'.format(master_seed))
  ins_loc_rng, ins_len_rng, base_sel_rng, het_rng, copy_rng = util.initialize_rngs(master_seed, 5)

  for chrom in chromosome:
    ref_seq = ref[chrom]
    ins_locs, = numpy.nonzero(ins_loc_rng.rand(len(ref_seq)) < p)
    if ins_locs.size == 0: continue  # No variants here
    ins_lens = ins_len_rng.randint(low=ins_len_lo, high=ins_len_hi+1, size=ins_locs.size)
    het_type = util.het(ins_locs.size, phet, het_rng, copy_rng)

    yield {chrom:[Variation(pos + 1, pos + 2, ref_seq[pos], ref_seq[pos] +
           base_sel_rng.choice(['A','C','G','T'], size=ins_len, replace=True, p=[.3, .2, .2, .3]).tostring(), het)
           for het, pos, ins_len in zip(het_type, ins_locs, ins_lens) if ref_seq[pos] != 'N']}

if __name__ == "__main__":
    print _description
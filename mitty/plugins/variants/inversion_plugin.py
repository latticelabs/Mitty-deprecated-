"""This is the stock inversion generator."""
from mitty.lib.util import initialize_rngs

__explain__ = """
Example parameter snippet:

    {
        "chromosome": [1],
        "model": "inversion",
        "phet": 0.5,
        "p": 0.01,
        "inv_len_lo": 10,
        "inv_len_hi": 100,
        "inv_loc_rng_seed": 1,
        "inv_len_rng_seed": 2,
        "het_rng_seed": 3,
        "copy_rng_seed": 4
    }
"""
import numpy
import logging
import string
import mitty.plugins.variants.util as util
from mitty.lib.variation import new_variation

logger = logging.getLogger(__name__)
DNA_complement = string.maketrans('ATCGNatcg', 'TAGCNtagc')

def variant_generator(ref={},
             chromosome=None,
             p=0.0001,
             phet=0.5,
             inv_len_lo=10,
             inv_len_hi=100,
             master_seed=None,
             inv_loc_rng_seed=1,
             inv_len_rng_seed=2,
             het_rng_seed=3,
             copy_rng_seed=4,
             **kwargs):

  try:
    vg = variant_generator
    inv_loc_rng, inv_len_rng, het_rng, copy_rng = vg.inv_loc_rng, vg.inv_len_rng, vg.het_rng, vg.copy_rng
    logger.debug('Using previous RNG states')
  except AttributeError:
    if master_seed is not None:
      inv_loc_rng_seed, inv_len_rng_seed, het_rng_seed, copy_rng_seed = \
      numpy.random.RandomState(seed=master_seed).randint(100000000, size=4)
      logger.debug('Used master seed to generate seeds {:d}, {:d}, {:d}, {:d}'.
                 format(inv_loc_rng_seed, inv_len_rng_seed, het_rng_seed, copy_rng_seed))

    inv_loc_rng, inv_len_rng, het_rng, copy_rng = initialize_rngs(inv_loc_rng_seed, inv_len_rng_seed, het_rng_seed, copy_rng_seed)
    vg = variant_generator
    vg.inv_loc_rng, vg.inv_len_rng, vg.het_rng, vg.copy_rng = inv_loc_rng, inv_len_rng, het_rng, copy_rng

  for chrom in chromosome:
    ref_seq = ref[chrom]
    inv_locs, = numpy.nonzero(inv_loc_rng.rand(len(ref_seq)) < p)
    inv_lens = inv_len_rng.randint(low=inv_len_lo, high=inv_len_hi+1, size=inv_locs.size)
    het_type = util.het(inv_locs.size, phet, het_rng, copy_rng)


    yield {chrom:[new_variation(pos + 1, pos + 1+inv_len, ref_seq[pos:pos + inv_len],
                ref_seq[pos:pos + inv_len].translate(DNA_complement)[::-1], het)
              for het, pos, inv_len in zip(het_type, inv_locs, inv_lens)
              if ('N'  not in ref_seq[pos:pos + inv_len]) and ('\n' not in ref_seq[pos:pos + inv_len])
                                                            and (pos + inv_len < len(ref_seq))]}

if __name__ == "__main__":
  print __explain__
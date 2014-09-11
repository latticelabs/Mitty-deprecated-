"""This is the stock insertion generator.

Note: This never generates a deletion at the first base of a sequence.

"""
__explain__ = """
Example parameter snippet:

    {
        "chromosome": [1],
        "phet": 0.5,
        "p": 0.01,
        "ins_len_lo": 100,
        "ins_len_hi": 10000,
        "ins_loc_rng_seed": 1,
        "ins_len_rng_seed": 2,
        "base_sel_rng_seed": 3,
        "het_rng_seed": 4,
        "copy_rng_seed": 5
    }
"""
import string
import numpy
import logging
import mitty.plugins.variants.util as util
from mitty.lib.variation import Variation


logger = logging.getLogger(__name__)


def variant_generator(ref={},
                      chromosome=None,
                      p=0.01,
                      phet=0.5,
                      ins_len_lo=100,
                      ins_len_hi=10000,
                      master_seed=None,
                      ins_loc_rng_seed=1,
                      ins_len_rng_seed=2,
                      base_sel_rng_seed=3,
                      het_rng_seed=4,
                      copy_rng_seed=5,
                      **kwargs):
  
  try:
      vg = variant_generator
      
      ins_loc_rng, ins_len_rng,base_sel_rng, het_rng, copy_rng = vg.ins_loc_rng, vg.ins_len_rng,vg.base_sel_rng, vg.het_rng, vg.copy_rng
      logger.debug('Using previous RNG states')
  
  except AttributeError:
      if master_seed is not None:

        ins_loc_rng_seed, ins_len_rng_seed, base_sel_rng_seed, het_rng_seed, copy_rng_seed = \
        numpy.random.RandomState(seed=master_seed).randint(100000000, size=5)
        logger.debug('Used master seed to generate seeds {:d}, {:d}, {:d}, {:d}, {:d}'.
                 format(ins_loc_rng_seed, ins_len_rng_seed, base_sel_rng_seed, het_rng_seed, copy_rng_seed))


      ins_loc_rng, ins_len_rng,base_sel_rng, het_rng, copy_rng= util.initialize_rngs(ins_loc_rng_seed, ins_len_rng_seed,
                                                                        base_sel_rng_seed, het_rng_seed, copy_rng_seed)
      vg = variant_generator
      vg.ins_loc_rng, vg.ins_len_rng,vg.base_sel_rng, vg.het_rng, vg.copy_rng= \
      ins_loc_rng, ins_len_rng,base_sel_rng, het_rng, copy_rng
      




  for chrom in chromosome:
    ref_seq = ref[chrom]
    ins_locs, = numpy.nonzero(ins_loc_rng.rand(len(ref_seq)) < p)
    if ins_locs.size == 0: continue  # No variants here
    ins_lens = ins_len_rng.randint(low=ins_len_lo, high=ins_len_hi+1, size=ins_locs.size)
    het_type = util.het(ins_locs.size, phet, het_rng, copy_rng)
    

    yield {chrom:[Variation(pos + 1, pos + 2, ref_seq[pos], ref_seq[pos] +
            base_sel_rng.choice(['A','C','G','T'], size=ins_len, replace=True, p=[.3, .2, .2, .3]).tostring(), het)
              for het, pos, ins_len in zip(het_type, ins_locs, ins_lens) if ref_seq[pos] != 'N' and ref_seq[pos] != '\n']}

if __name__ == "__main__":
    print __explain__


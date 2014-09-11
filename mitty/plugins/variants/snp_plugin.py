"""This is the stock SNP plugin. It uses four independent RNGs to locate SNPs along a reference sequence, assign each
SNP a zygosity and assign an ALT base.
"""
import numpy
import mitty.plugins.variants.util as util
from mitty.lib.variation import Variation
import logging
logger = logging.getLogger(__name__)

__example_param_text = """
{
  "chromosome": [1],      # List of chromosomes to apply the variant to
  "p": 0.001,             # probability that the SNP will happen at any given base
  "phet": 0.5,            # probability that the variant will be heterozygous
  "base_loc_rng_seed": 1, # Seeds for the RNGs
  "base_sub_rng_seed": 2,
  "het_rng_seed": 3,
  "copy_rng_seed": 4
}
"""

_description = """
This is the stock SNP plugin. A typical parameter set resembles
""" + __example_param_text

_example_params = eval(__example_param_text)


# This is a substitution base substitution table
base_sub_dict = {
  'A': ['T', 'C', 'G'],  # 'A' -> 'TCG'
  'C': ['G', 'A', 'T'],  # 'C' -> 'GAT'
  'G': ['A', 'T', 'C'],  # 'G' -> 'ATC'
  'T': ['C', 'G', 'A'],  # 'T' -> 'CGA'
  'N': ['N', 'N', 'N']

}


def variant_generator(ref={},
                      chromosome=None,
                      phet=0.5,
                      p=0.01,
                      master_seed=None,
                      base_loc_rng_seed=1,
                      base_sub_rng_seed=2,
                      het_rng_seed=3,
                      copy_rng_seed=4,
                      **kwargs):
  try:
    vg = variant_generator
    
    base_loc_rng, base_sub_rng, het_rng, copy_rng = vg.base_loc_rng, vg.base_sub_rng, vg.het_rng, vg.copy_rng
    logger.debug('Using previous RNG states')
  except AttributeError:
    if master_seed is not None:
      base_loc_rng_seed, base_sub_rng_seed, het_rng_seed, copy_rng_seed = \
        numpy.random.RandomState(seed=master_seed).randint(100000000, size=4)
      logger.debug('Used master seed to generate seeds {:d}, {:d}, {:d}, {:d}'.
                   format(base_loc_rng_seed, base_sub_rng_seed, het_rng_seed, copy_rng_seed))

    base_loc_rng, base_sub_rng, het_rng, copy_rng = util.initialize_rngs(base_loc_rng_seed, base_sub_rng_seed, het_rng_seed, copy_rng_seed)
    vg = variant_generator
    vg.base_loc_rng, vg.base_sub_rng, vg.het_rng, vg.copy_rng = base_loc_rng, base_sub_rng, het_rng, copy_rng

  for chrom in chromosome:
    ref_chrom = ref[chrom]
    snp_locs = util.place_poisson(base_loc_rng, p, len(ref_chrom))
    het_type = util.het(snp_locs.size, phet, het_rng, copy_rng)
    base_subs = base_sub_rng.randint(3, size=snp_locs.size)


    yield {chrom: [Variation(pos + 1, pos + 2, ref_chrom[pos], base_sub_dict[ref_chrom[pos].upper()][bs], het)
                   for pos, bs, het in zip(snp_locs, base_subs, het_type) if ref_chrom[pos] != 'N' and ref_chrom[pos] != '\n']}
    # +1 because VCF files are 1 indexed
    #alts will be 0 if ref is not one of ACTG


def test():
  ref = {
    1: 'ACTGACTGACTG',
    2: 'ACTGACTGACTGACTGACTGACTG'
  }
  vg = variant_generator(ref, chromosome=[1, 2], p=1.0)
  for v_list in vg:
    assert isinstance(v_list.values()[0][0], Variation), v_list

if __name__ == "__main__":
  print _description
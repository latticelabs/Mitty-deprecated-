"""This is the stock SNP plugin. It uses four independent RNGs to locate SNPs along a reference sequence, assign each
SNP a zygosity and assign an ALT base.
"""
import numpy
import mitty.Plugins.variants.util as util
from mitty.variation import Variation
import logging
logger = logging.getLogger(__name__)

__explain__ = """
Example parameter snippet:

    {
        "chromosome": [1],
        "model": "snp",
        "phet": 0.5,
        "p": 0.01,
        "master_seed": 0,
        "base_loc_rng_seed": 1,
        "base_sub_rng_seed": 2,
        "het_rng_seed": 3,
        "copy_rng_seed": 4
    }

Notes:
1. If you add a master_seed the other seeds are not used. Supplying a master seed to mutate will override this


In order to work with mutate.py the model needs a method called variants that returns three lists

description, footprint, vcf_line

higher order description  ()  - as needed to describe complex variants - will make VCF more sophisticated as needed later)
footprints                (het, chrom, pos_st, pos_nd) ... used for collision detection
vcf list                  (chrom, pos, id, ref, alt, qual, filter, info, format, sample) ... (as many as needed) - for simple VCF file


1 -> variant on copy 1
2 -> variant on copy 2
3 -> variant on both copies

"""
# This is a substitution base substitution table
base_sub_mat = numpy.zeros((256, 3), dtype='u1')
base_sub_mat[65, :] = [ord(c) for c in ['T', 'C', 'G']]  # 'A' -> 'TCG'
base_sub_mat[67, :] = [ord(c) for c in ['G', 'A', 'T']]  # 'C' -> 'GAT'
base_sub_mat[71, :] = [ord(c) for c in ['A', 'T', 'C']]  # 'G' -> 'ATC'
base_sub_mat[84, :] = [ord(c) for c in ['C', 'G', 'A']]  # 'T' -> 'CGA'


def variant_generator(ref_fp=None,
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
    ref_seq = ref_fp['sequence/{:d}/1'.format(chrom)][:]
    #snp_locs, = numpy.nonzero(base_loc_rng.rand(len(ref_seq)) < p)
    snp_locs = util.place_poisson(base_loc_rng, p, ref_seq.size)
    het_type = util.het(snp_locs.size, phet, het_rng, copy_rng)
    base_subs = base_sub_rng.randint(3, size=snp_locs.size)
    if snp_locs.size:
      refs = ref_seq[snp_locs]
      refs_s = refs.tostring()
      alts = base_sub_mat[refs, base_subs]
      alts_s = alts.tostring()
    else:
      refs, refs_s, alts, alts_s = [], [], [], []

    yield {chrom: [Variation(pos + 1, pos + 2, ref, alt, het) for pos, ref, alt, het, valid in zip(snp_locs, refs_s, alts_s, het_type, alts) if valid]}
    # +1 because VCF files are 1 indexed
    #alts will be 0 if ref is not one of ACTG

if __name__ == "__main__":
  print __explain__
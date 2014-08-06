"""This is the stock SNP plugin. It uses four independent RNGs to locate SNPs along a reference sequence, assign each
SNP a zygosity and assign an ALT base.
"""
import numpy
import logging
import mitty.Plugins.variants.util as util
from mitty.variation import Variation

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

logger = logging.getLogger(__name__)
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
    ref_seq = ref_fp['sequence/{:d}/1'.format(chrom)][:]#.tostring()  # Very cheap operation
    snp_locs, = numpy.nonzero(base_loc_rng.rand(len(ref_seq)) < p)
    het_type = util.het(snp_locs.size, phet, het_rng, copy_rng)
    base_subs = base_sub_rng.randint(3, size=snp_locs.size)

    refs = ref_seq[snp_locs]
    refs_s = refs.tostring()
    alts = base_sub_mat[refs, base_subs]
    alts_s = alts.tostring()

    yield {chrom: [Variation(pos + 1, pos + 2, ref, alt, het) for pos, ref, alt, het, valid in zip(snp_locs, refs_s, alts_s, het_type, alts) if valid]}
    # +1 because VCF files are 1 indexed

def variants(ref_fp=None,
             chromosome=None,
             phet=0.5,
             p=0.01,
             master_seed=None,
             base_loc_rng_seed=1,
             base_sub_rng_seed=2,
             het_rng_seed=3,
             copy_rng_seed=4,
             **kwargs):
  """This is the stock SNP generator and returns snp locations in
  a poisson distributed fashion.
  Inputs:
    ref_seq              - The reference sequence
    ref_seq_len          - length of whole sequence (needed to compute start and stop)
    phet                 - probability of having heterozygous mutation
    p                    - probability of SNPs
    start_snps_frac      - start generating snps from here (0.0, 1.0)
    stop_snps_frac       - stop generating snps after this (0.0, 1.0) stop_snps_frac > start_snps_frac
    het_rng_seed         - rng used to decide if genotype is heterozygous or not (0/1 or 1/0  vs 1/1)
    strand_rng_seed      - rng used to decide which het type 0/1 or 1/0
    poisson_rng_seed     - SNP locator rng numpy.random.RandomState(seed)
    base_sub_rng_seed    - rng used to select ALT bases
    kwargs               - absorbs any other parameters it does not use

  Outputs:
    variant              - (POS, REF, ALT, GT, skip, list(footprints))

  Any section with anything other than GATC should be skipped
  >>> args = {'p': .1, 'poisson_rng_seed': 1, 'base_sub_rng_seed': 2}; \
  ref_seq='ACGTACGTAcGTACGNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'; ref_seq_len = len(ref_seq); \
  gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=100, **args);
  >>> for n in range(10): print next(gen,None)
  (22, 'G', 'A', '1/1', 24, None)
  (31, 'T', 'A', '1/1', 33, None)
  (40, 'A', 'G', '1/1', 42, None)
  (49, 'C', 'G', '1/1', 51, None)
  (55, 'T', 'A', '1/1', 57, None)
  None
  None
  None
  None
  None

  Test with one block
  >>> args = {'p': .1, 'poisson_rng_seed': 1, 'base_sub_rng_seed': 2}; \
  ref_seq='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'; ref_seq_len = len(ref_seq); \
  gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=100, **args);
  >>> for n in range(10): print next(gen,None)
  (9, 'C', 'G', '1/1', 11, None)
  (15, 'T', 'G', '1/1', 17, None)
  (22, 'G', 'A', '1/1', 24, None)
  (31, 'T', 'A', '1/1', 33, None)
  (40, 'A', 'G', '1/1', 42, None)
  (49, 'C', 'G', '1/1', 51, None)
  (55, 'T', 'A', '1/1', 57, None)
  None
  None
  None

  Test with multiple blocks - should be same answer even though we have changed the block size
  >>> gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=1, **args);
  >>> for n in range(10): print next(gen,None)
  (9, 'C', 'G', '1/1', 11, None)
  (15, 'T', 'G', '1/1', 17, None)
  (22, 'G', 'A', '1/1', 24, None)
  (31, 'T', 'A', '1/1', 33, None)
  (40, 'A', 'G', '1/1', 42, None)
  (49, 'C', 'G', '1/1', 51, None)
  (55, 'T', 'A', '1/1', 57, None)
  None
  None
  None

  Test with heterozygosity
  >>> args = {'phet': 0.5, 'p': .1, 'het_rng_seed': 3, 'strand_rng_seed': 5, 'poisson_rng_seed': 1, 'base_sub_rng_seed': 2}; \
  ref_seq='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'; ref_seq_len = len(ref_seq); \
  gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=100, **args);
  >>> for n in range(10): print next(gen,None)
  (9, 'C', 'G', '1/1', 11, None)
  (15, 'T', 'G', '1/1', 17, None)
  (22, 'G', 'A', '1/0', 24, None)
  (31, 'T', 'A', '1/1', 33, None)
  (40, 'A', 'G', '1/1', 42, None)
  (49, 'C', 'G', '1/1', 51, None)
  (55, 'T', 'A', '0/1', 57, None)
  None
  None
  None
  """
  if master_seed:
    base_loc_rng_seed, base_sub_rng_seed, het_rng_seed, copy_rng_seed = \
      numpy.random.RandomState(seed=master_seed).randint(100000000, size=4)
    logger.debug('Used master seed to generate seeds {:d}, {:d}, {:d}, {:d}'.
                 format(base_loc_rng_seed, base_sub_rng_seed, het_rng_seed, copy_rng_seed))

  base_loc_rng = numpy.random.RandomState(seed=base_loc_rng_seed)
  base_sub_rng = numpy.random.RandomState(seed=base_sub_rng_seed)
  het_rng = numpy.random.RandomState(seed=het_rng_seed)
  copy_rng = numpy.random.RandomState(seed=copy_rng_seed)

  description, footprint, vcf_line = [], [], []

  for chrom in chromosome:
    ref_seq = ref_fp['sequence/{:d}/1'.format(chrom)][:].tostring()  # Very cheap operation
    base_locs, = numpy.nonzero(base_loc_rng.rand(len(ref_seq)) < p)
    het_type = numpy.empty((base_locs.size,), dtype='u1')
    het_type.fill(3)  # Homozygous
    idx_het, = numpy.nonzero(het_rng.rand(het_type.size) < phet)  # Heterozygous locii
    het_type[idx_het] = 1  # On copy 1
    het_type[idx_het[numpy.nonzero(copy_rng.rand(idx_het.size) < 0.5)[0]]] = 2  # On copy 2
    base_subs = base_sub_rng.randint(3, size=base_locs.size)

    for het, pos, bsub in zip(het_type, base_locs, base_subs):
      ref = ref_seq[pos]
      if ref in ['G', 'A', 'T', 'C']:
        alt = base_sub_mat[ref][bsub]
        gt = gt_string[het]
      else:
        continue  # Not valid, skip to next

      description.append('SNP')
      footprint.append([(het, chrom-1, pos, pos + 1)])  # Chrom is internal numbering, starts from 0
      # [(het, chrom, pos_st, pos_nd)]
      vcf_line.append([(chrom, pos+1, '.', ref, alt, 100, 'PASS', '.', 'GT', gt)])  # POS is VCF number starts from 1
      # CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample
  return description, footprint, vcf_line

if __name__ == "__main__":
  print __explain__
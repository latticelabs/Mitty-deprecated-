"""This is the stock insertion generator.

Note: This never generates a deletion at the first base of a sequence.

"""
__explain__ = """
Example parameter snippet:

    {
        "chromosome": [1],
        "model": "insert",
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

logger = logging.getLogger(__name__)
#             0      1      2      3
gt_string = ['0/0', '0/1', '1/0', '1/1']  # The types of genotypes


def variants(ref_fp=None,
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
  """A generator which returns a variant when asked for. This is the stock SNP generator and returns snp locations in
  a poisson distributed fashion.
  Inputs:
    ref_seq              - The reference sequence
    ref_seq_len          - length of whole sequence (needed to compute start and stop)
    p                    - probability of inserts
    phet                 - probability of having heterozygous mutation
    dist                 - The kind of distribution for inserts
                 {'type': 'uniform', 'low': 100, 'high': 10000} generate insertions uniformly between these lengths
                 {'type': 'geometric', 'a': 0.1}
    copy_rng_seed        - rng used to decide which het type 0/1 or 1/0
    ins_loc_rng_seed     - INS locator rng -> numpy.random.RandomState(seed)
    ins_len_rng_seed     - rng used to determine length of insert
    base_sel_rng_seed    - seed for rng used to determine bases to insert
    het_rng_seed         - rng used to decide if genotype is heterozygous or not (0/1 or 1/0  vs 1/1)
    kwargs               - absorbs any other parameters it does not use

  Outputs:
    variant              - description, footprint, vcf_line


  Test with 'N' regions. There should be no insertions after an 'N'
  >>> args = {'p_ins': .1, 'lam_ins': 3, 'ins_loc_rng_seed': 1, 'ins_len_rng_seed': 2, 'base_sel_rng_seed': 3}; \
  ref_seq='ACGTACGTANGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'; ref_seq_len = len(ref_seq); \
  gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=100, **args);
  >>> for n in range(10): print next(gen,None)
  (15, 'T', 'TGAC', '1/1', 17, None)
  (22, 'G', 'GGG', '1/1', 24, None)
  (31, 'T', 'TGAA', '1/1', 33, None)
  (40, 'A', 'ACTCA', '1/1', 42, None)
  (49, 'C', 'CAT', '1/1', 51, None)
  (55, 'T', 'TA', '1/1', 57, None)
  None
  None
  None
  None


  Test with one block
  >>> args = {'p_ins': .1, 'lam_ins': 3, 'ins_loc_rng_seed': 1, 'ins_len_rng_seed': 2, 'base_sel_rng_seed': 3}; \
  ref_seq='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'; ref_seq_len = len(ref_seq); \
  gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=100, **args);
  >>> for n in range(10): print next(gen,None)
  (9, 'C', 'CT', '1/1', 11, None)
  (15, 'T', 'TGAC', '1/1', 17, None)
  (22, 'G', 'GGG', '1/1', 24, None)
  (31, 'T', 'TGAA', '1/1', 33, None)
  (40, 'A', 'ACTCA', '1/1', 42, None)
  (49, 'C', 'CAT', '1/1', 51, None)
  (55, 'T', 'TA', '1/1', 57, None)
  None
  None
  None

  Test with multiple blocks - should be same answer even though we have changed the block size
  >>> gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=1, **args);
  >>> for n in range(10): print next(gen,None)
  (9, 'C', 'CT', '1/1', 11, None)
  (15, 'T', 'TGAC', '1/1', 17, None)
  (22, 'G', 'GGG', '1/1', 24, None)
  (31, 'T', 'TGAA', '1/1', 33, None)
  (40, 'A', 'ACTCA', '1/1', 42, None)
  (49, 'C', 'CAT', '1/1', 51, None)
  (55, 'T', 'TA', '1/1', 57, None)
  None
  None
  None

  Test with heterozygosity
  >>> args = {'phet': 0.5, 'p_ins': .1, 'lam_ins': 3, 'het_rng_seed': 3, 'strand_rng_seed': 5, 'ins_loc_rng_seed': 1, 'ins_len_rng_seed': 2, 'base_sel_rng_seed': 3}; \
  ref_seq='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'; ref_seq_len = len(ref_seq); \
  gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=100, **args);
  >>> for n in range(10): print next(gen,None)
  (9, 'C', 'CT', '1/1', 11, None)
  (15, 'T', 'TGAC', '1/1', 17, None)
  (22, 'G', 'GGG', '1/0', 24, None)
  (31, 'T', 'TGAA', '1/1', 33, None)
  (40, 'A', 'ACTCA', '1/1', 42, None)
  (49, 'C', 'CAT', '1/1', 51, None)
  (55, 'T', 'TA', '0/1', 57, None)
  None
  None
  None
  """
  if master_seed:
    ins_loc_rng_seed, ins_len_rng_seed, base_sel_rng_seed, het_rng_seed, copy_rng_seed = \
      numpy.random.RandomState(seed=master_seed).randint(100000000, size=5)
    logger.debug('Used master seed to generate seeds {:d}, {:d}, {:d}, {:d}, {:d}'.
                 format(ins_loc_rng_seed, ins_len_rng_seed, base_sel_rng_seed, het_rng_seed, copy_rng_seed))

  ins_loc_rng = numpy.random.RandomState(seed=ins_loc_rng_seed)
  ins_len_rng = numpy.random.RandomState(seed=ins_len_rng_seed)
  base_sel_rng = numpy.random.RandomState(seed=base_sel_rng_seed)
  het_rng = numpy.random.RandomState(seed=het_rng_seed)
  copy_rng = numpy.random.RandomState(seed=copy_rng_seed)

  description, footprint, vcf_line = [], [], []

  for chrom in chromosome:
    ref_seq = ref_fp['sequence/{:d}/1'.format(chrom)][:].tostring()  # Very cheap operation
    ins_locs, = numpy.nonzero(ins_loc_rng.rand(len(ref_seq)) < p)
    if ins_locs.size == 0: continue  # No variants here
    ins_lens = ins_len_rng.randint(low=ins_len_lo, high=ins_len_hi+1, size=ins_locs.size)
    het_type = numpy.empty((ins_locs.size,), dtype='u1')
    het_type.fill(3)  # Homozygous
    idx_het, = numpy.nonzero(het_rng.rand(het_type.size) < phet)  # Heterozygous locii
    het_type[idx_het] = 1  # On copy 1
    het_type[idx_het[numpy.nonzero(copy_rng.rand(idx_het.size) < 0.5)[0]]] = 2  # On copy 2
    for het, pos, ins_len in zip(het_type, ins_locs, ins_lens):
      ref = ref_seq[pos]
      if ref == 'N':
        continue  # Not valid, skip to next
      else:
        alt = ref + base_sel_rng.choice(['A','C','G','T'], size=ins_len, replace=True, p=[.3, .2, .2, .3]).tostring()
        gt = gt_string[het]

      description.append('Insert')
      footprint.append([(het, chrom-1, pos, pos + 1)])  # Chrom is internal numbering, starts from 0
      # [(het, chrom, pos_st, pos_nd)]
      vcf_line.append([(chrom, pos+1, '.', ref, alt, 100, 'PASS', '.', 'GT', gt)])  # POS is VCF number starts from 1
      # CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample

    logger.debug('Generated {:d} insertions with min={:d},mean={:f},max={:d}'.
                 format(len(vcf_line), int(ins_lens.min()), ins_lens.mean(), int(ins_lens.max())))
  return description, footprint, vcf_line
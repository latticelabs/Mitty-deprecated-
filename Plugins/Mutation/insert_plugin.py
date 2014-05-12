"""This is the stock insertion generator. Please see Readme for details on how to write variant generator plugins for
mutate.py.

The plugin uses three random number generators. The first creates poisson distributed numbers to locate the deletions.
The second creates poisson distributed numbers to determine deletion lengths. This is equivalent to modeling deletion
termination as a bernoulli process. The third is a uniform random number generator that is used to generate bases for
the insertions. The plugin generates all bases at the same rate.

Note: This never generates a deletion at the first base of a sequence."""
import numpy
import logging

logger = logging.getLogger(__name__)
bases = 'GATC'
het_type = ['0/1', '1/0']  # The two types of het

def variant(ref_seq=None, ref_seq_len=0,
            phet=0.0, p_ins=0.01, lam_ins=5,
            start_ins_frac=0.0, stop_ins_frac=1.0,
            het_rng_seed=1,
            strand_rng_seed=4,
            ins_loc_rng_seed=0,
            ins_len_rng_seed=1,
            base_sel_rng_seed=2,
            block_size=10000, **kwargs):
  """A generator which returns a variant when asked for. This is the stock SNP generator and returns snp locations in
  a poisson distributed fashion.
  Inputs:
    ref_seq              - The reference sequence
    ref_seq_len          - length of whole sequence (needed to compute start and stop)
    phet                 - probability of having heterozygous mutation
    p_ins                - probability of inserts
    lam_ins              - mean length of poisson distributed insert lengths
    start_ins_frac       - start generating inserts from here (0.0, 1.0)
    stop_ins_frac        - stop generating inserts after this (0.0, 1.0) stop_ins_frac > start_ins_frac
    het_rng_seed         - rng used to decide if genotype is heterozygous or not (0/1 or 1/0  vs 1/1)
    strand_rng_seed      - rng used to decide which het type 0/1 or 1/0
    ins_loc_rng_seed     - INS locator rng -> numpy.random.RandomState(seed)
    ins_len_rng_seed     - rng used to determine length of insert
    base_sel_rng_seed    - seed for rng used to determine bases to insert
    block_size           - how many random numbers should we geenrate at a time
    kwargs               - absorbs any other parameters it does not use

  Outputs:
    variant              - (POS, REF, ALT, GT, skip, list(footprints))


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
  def get_locs_and_lens():  # Simply a convenience.
    het_or_not = het_rng.rand(block_size)
    strand_no = strand_rng.randint(2, size=block_size)
    loc_diff = numpy.maximum(ins_loc_rng.poisson(lam=1.0 / p_ins, size=block_size), 1)  # Please sir, can I have some more?
    ins_lens = ins_len_rng.poisson(lam=lam_ins, size=block_size)
    bc = base_sel_rng.randint(0, 4, size=numpy.sum(ins_lens))
    csl = numpy.concatenate(([0],numpy.cumsum(ins_lens)))
    base_calls = [bc[csl[n]:csl[n+1]] for n in range(block_size)]
    locs = numpy.cumsum(loc_diff) + start_offset
    return het_or_not, strand_no, locs, ins_lens, locs[-1], base_calls

  het_rng = numpy.random.RandomState(seed=het_rng_seed)
  strand_rng = numpy.random.RandomState(seed=strand_rng_seed)
  ins_loc_rng = numpy.random.RandomState(seed=ins_loc_rng_seed)
  ins_len_rng = numpy.random.RandomState(seed=ins_len_rng_seed)
  base_sel_rng = numpy.random.RandomState(seed=base_sel_rng_seed)

  start_offset = int(ref_seq_len * start_ins_frac)
  ins_end = int(ref_seq_len * stop_ins_frac)

  het_or_not, strand_no, locs, ins_lens, start_offset, base_calls = get_locs_and_lens()
  internal_cntr = 0
  while locs[internal_cntr] < ins_end:
    vl = locs[internal_cntr]
    ref = ref_seq[vl]
    if ref == 'N':
      alt = None  # Don't do insertions is regions with Ns
    else:
      alt = ref + ''.join([bases[k] for k in base_calls[internal_cntr]])
      gt = '1/1' if het_or_not[internal_cntr] > phet else het_type[strand_no[internal_cntr]]

    internal_cntr += 1
    if internal_cntr == locs.size:
      het_or_not, strand_no, locs, ins_lens, start_offset, base_calls = get_locs_and_lens()
      internal_cntr = 0

    if alt is not None:
      yield (vl, ref, alt, gt, vl+2, None)  # POS, REF, ALT, skipto, list(footprints)
                                         # footprints, in this case, is None, since we simply skip forward
      # We have vl + 2 because we want 1 base buffer between variants (see Readme)


if __name__ == "__main__":
  import doctest
  doctest.testmod()

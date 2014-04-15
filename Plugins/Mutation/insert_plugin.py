"""TODO: Write docs"""
import numpy
import logging

logger = logging.getLogger(__name__)
bases = 'GATC'
# base_occurrence_prob = {  # How likely are we to get a base? Used for generating insertions
#                   'G': .2,
#                   'A': .3,
#                   'T': .3,
#                   'C': .2
# }


def variant(p_ins=0.01, lam_ins=5,
            ref_seq=None, ref_seq_len=0, start_ins_frac=0.0, stop_ins_frac=1.0,
            ins_loc_rng_seed=0,
            ins_len_rng_seed=1,
            base_sel_rng_seed=2,
            block_size=10000, **kwargs):
  """A generator which returns a variant when asked for. This is the stock SNP generator and returns snp locations in
  a poisson distributed fashion.
  Inputs:
    ref_seq              - The reference sequence
    ref_seq_len          - length of whole sequence (needed to compute start and stop)
    p_ins                - probability of inserts
    lam_ins              - mean length of poisson distributed insert lengths
    start_ins_frac       - start generating inserts from here (0.0, 1.0)
    stop_ins_frac        - stop generating inserts after this (0.0, 1.0) stop_ins_frac > start_ins_frac
    ins_loc_rng_seed     - INS locator rng -> numpy.random.RandomState(seed)
    ins_len_rng_seed     - rng used to determine length of insert
    base_sel_rng_seed    - seed for rng used to determine bases to insert
    block_size           - how many random numbers should we geenrate at a time
    kwargs               - absorbs any other parameters it does not use

  Outputs:
    variant              - (POS, REF, ALT, skip, list(footprints))


  Test with one block
  >>> args = {'p_del': .1, 'lam_del': 3, 'del_loc_rng_seed': 1, 'del_len_rng_seed': 2}; \
  ref_seq='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'; ref_seq_len = len(ref_seq); \
  gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=100, **args);
  >>> for n in range(10): print next(gen,None)
  (9, 'CG', 'C', 11, None)
  (15, 'TACG', 'T', 19, None)
  (22, 'GTA', 'G', 25, None)
  (31, 'TACG', 'T', 35, None)
  (40, 'ACGTA', 'A', 45, None)
  (49, 'CGT', 'C', 52, None)
  (55, 'TA', 'T', 57, None)
  None
  None
  None

  Test with multiple blocks - should be same answer even though we have changed the block size
  >>> gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=1, **args);
  >>> for n in range(10): print next(gen,None)
  (9, 'CG', 'C', 11, None)
  (15, 'TACG', 'T', 19, None)
  (22, 'GTA', 'G', 25, None)
  (31, 'TACG', 'T', 35, None)
  (40, 'ACGTA', 'A', 45, None)
  (49, 'CGT', 'C', 52, None)
  (55, 'TA', 'T', 57, None)
  None
  None
  None
  """
  def get_locs_and_lens():  # Simply a convenience.
    loc_diff = numpy.maximum(ins_loc_rng.poisson(lam=1.0 / p_ins, size=block_size), 1)  # Please sir, can I have some more?
    ins_lens = ins_len_rng.poisson(lam=lam_ins, size=block_size)
    bc = base_sel_rng.randint(0, 4, size=numpy.sum(ins_lens))
    csl = numpy.concatenate(([0],numpy.cumsum(ins_lens)))
    base_calls = [bc[csl[n]:csl[n+1]] for n in range(block_size)]
    locs = numpy.cumsum(loc_diff) + start_offset
    return locs, ins_lens, locs[-1], base_calls

  ins_loc_rng = numpy.random.RandomState(seed=ins_loc_rng_seed)
  ins_len_rng = numpy.random.RandomState(seed=ins_len_rng_seed)
  base_sel_rng = numpy.random.RandomState(seed=base_sel_rng_seed)

  start_offset = int(ref_seq_len * start_ins_frac)
  ins_end = int(ref_seq_len * stop_ins_frac)

  locs, ins_lens, start_offset, base_calls = get_locs_and_lens()
  internal_cntr = 0
  while locs[internal_cntr] < ins_end:
    vl = locs[internal_cntr]
    ref = ref_seq[vl]
    if 'N' in ref:
      alt = None  # Don't do insertions is regions with Ns
    else:
      alt = ref + ''.join([bases[k] for k in base_calls[internal_cntr]])

    internal_cntr += 1
    if internal_cntr == locs.size:
      locs, ins_lens, start_offset, base_calls = get_locs_and_lens()
      internal_cntr = 0

    if alt is not None:
      yield (vl, ref, alt, vl+1, None)  # POS, REF, ALT, skipto, list(footprints)
                                         # footprints, in this case, is None, since we simply skip forward


if __name__ == "__main__":
  import doctest
  doctest.testmod()

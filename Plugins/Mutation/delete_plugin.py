"""TODO: Write docs"""
import numpy
import logging

logger = logging.getLogger(__name__)


def variant(p_del=0.01, lam_del=5,
            ref_seq=None, ref_seq_len=0, start_dels_frac=0.0, stop_dels_frac=1.0,
            del_loc_rng_seed=1,
            del_len_rng_seed=1,
            block_size=10000, **kwargs):
  """A generator which returns a variant when asked for. This is the stock SNP generator and returns snp locations in
  a poisson distributed fashion.
  Inputs:
    ref_seq              - The reference sequence
    ref_seq_len          - length of whole sequence (needed to compute start and stop)
    p_del                - probability of deletes
    lam_del              - mean length of poisson distributed delete lengths
    start_dels_frac      - start generating dels from here (0.0, 1.0)
    stop_dels_frac       - stop generating dels after this (0.0, 1.0) stop_snps_frac > start_snps_frac
    poisson_rng_seed     - SNP locator rng numpy.random.RandomState(seed)
    del_len_rng_seed     - rng used to determine length of delete
    kwargs               - absorbs any other parameters it does not use

  Outputs:
    variant              - (POS, REF, ALT, skip, list(footprints))


  Test with 'N's. No deletions should straddle a region with N
  >>> args = {'p_del': .1, 'lam_del': 3, 'del_loc_rng_seed': 1, 'del_len_rng_seed': 2}; \
  ref_seq='ACGTACGTANGTACGTACGTACGTACGTACGTACGTACGTACNTACGTACGTACGTACGT'; ref_seq_len = len(ref_seq); \
  gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=100, **args);
  >>> for n in range(10): print next(gen,None)
  (15, 'TACG', 'T', 19, None)
  (22, 'GTA', 'G', 25, None)
  (31, 'TACG', 'T', 35, None)
  (49, 'CGT', 'C', 52, None)
  (55, 'TA', 'T', 57, None)
  None
  None
  None
  None
  None


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
    loc_diff = numpy.maximum(del_loc_rng.poisson(lam=1.0 / p_del, size=block_size), 1)  # Please sir, can I have some more?
    del_lens = del_len_rng.poisson(lam=lam_del, size=block_size)
    locs = numpy.cumsum(loc_diff) + start_offset
    return locs, del_lens, locs[-1]

  del_loc_rng = numpy.random.RandomState(seed=del_loc_rng_seed)
  del_len_rng = numpy.random.RandomState(seed=del_len_rng_seed)

  start_offset = int(ref_seq_len * start_dels_frac)
  del_end = int(ref_seq_len * stop_dels_frac)

  locs, del_lens, start_offset = get_locs_and_lens()
  internal_cntr = 0
  while locs[internal_cntr] < del_end:
    vl = locs[internal_cntr]
    dl = del_lens[internal_cntr]
    ref = ref_seq[vl:vl+dl+1]
    if 'N' in ref:
      alt = None  # Very conservative - we only do deletions in completely known regions
    else:
      alt = ref[0]

    internal_cntr += 1
    if internal_cntr == locs.size:
      locs, del_lens, start_offset = get_locs_and_lens()
      internal_cntr = 0

    if alt is not None:
      yield (vl, ref, alt, vl + dl+1, None)  # POS, REF, ALT, skipto, list(footprints)
                                           # footprints, in this case, is None, since we simply skip forward


if __name__ == "__main__":
  import doctest
  doctest.testmod()

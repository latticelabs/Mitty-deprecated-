"""TODO: Write docs"""
import numpy
import logging

logger = logging.getLogger(__name__)
base_sub_mat = {  # GATC
                  'G': 'ATC',
                  'A': 'TCG',
                  'T': 'CGA',
                  'C': 'GAT'
}


def variant(p=0.01,
            ref_seq=None, ref_seq_len=0, start_snps_frac=0.0, stop_snps_frac=1.0,
            poisson_rng_seed=1,
            base_sub_rng_seed=1,
            block_size=10000, **kwargs):
  """A generator which returns a variant when asked for. This is the stock SNP generator and returns snp locations in
  a poisson distributed fashion.
  Inputs:
    ref_seq              - The reference sequence
    ref_seq_len          - length of whole sequence (needed to compute start and stop)
    p                    - probability of SNPs
    start_snps_frac      - start generating snps from here (0.0, 1.0)
    stop_snps_frac       - stop generating snps after this (0.0, 1.0) stop_snps_frac > start_snps_frac
    poisson_rng_seed     - SNP locator rng numpy.random.RandomState(seed)
    base_sub_rng_seed    - rng used to select ALT bases
    kwargs               - absorbs any other parameters it does not use

  Outputs:
    variant              - (POS, REF, ALT, skip, list(footprints))

  Any section with anything other than GATC should be skipped
  >>> args = {'p': .1, 'poisson_rng_seed': 1, 'base_sub_rng_seed': 2}; \
  ref_seq='ACGTACGTAcGTACGNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'; ref_seq_len = len(ref_seq); \
  gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=100, **args);
  >>> for n in range(10): print next(gen,None)
  (22, 'G', 'A', 23, None)
  (31, 'T', 'A', 32, None)
  (40, 'A', 'G', 41, None)
  (49, 'C', 'G', 50, None)
  (55, 'T', 'A', 56, None)
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
  (9, 'C', 'G', 10, None)
  (15, 'T', 'G', 16, None)
  (22, 'G', 'A', 23, None)
  (31, 'T', 'A', 32, None)
  (40, 'A', 'G', 41, None)
  (49, 'C', 'G', 50, None)
  (55, 'T', 'A', 56, None)
  None
  None
  None

  Test with multiple blocks - should be same answer even though we have changed the block size
  >>> gen = variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=1, **args);
  >>> for n in range(10): print next(gen,None)
  (9, 'C', 'G', 10, None)
  (15, 'T', 'G', 16, None)
  (22, 'G', 'A', 23, None)
  (31, 'T', 'A', 32, None)
  (40, 'A', 'G', 41, None)
  (49, 'C', 'G', 50, None)
  (55, 'T', 'A', 56, None)
  None
  None
  None

  """
  def get_locs_and_subs():  # Simply a convenience.
    poiss = numpy.maximum(poisson_rng.poisson(lam=1.0 / p, size=block_size), 1)  # Please sir, can I have some more?
    base_subs = base_sub_rng.randint(3, size=block_size)
    locs = numpy.cumsum(poiss) + start_offset
    return locs, base_subs, locs[-1]

  poisson_rng = numpy.random.RandomState(seed=poisson_rng_seed)
  base_sub_rng = numpy.random.RandomState(seed=base_sub_rng_seed)

  start_offset = int(ref_seq_len * start_snps_frac)
  snp_end = int(ref_seq_len * stop_snps_frac)

  locs, subs, start_offset = get_locs_and_subs()
  internal_cntr = 0
  while locs[internal_cntr] < snp_end:
    vl = locs[internal_cntr]
    ref = ref_seq[vl]
    if ref in ['G', 'A', 'T', 'C']:
      alt = base_sub_mat[ref][subs[internal_cntr]]
    else:
      alt = None  # Not valid, skip to next

    internal_cntr += 1
    if internal_cntr == locs.size:
      locs, subs, start_offset = get_locs_and_subs()
      internal_cntr = 0

    if alt is not None:
      yield (vl, ref, alt, vl + 2, None)  # POS, REF, ALT, skipto, list(footprints)
                                          # footprints, in this case, is None, since we simply skip forward
      # We have vl + 2 because we want 1 base buffer between variants, even SNPs (see Readme)

if __name__ == "__main__":
  import doctest
  doctest.testmod()

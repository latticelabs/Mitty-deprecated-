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


def candidate_variants(chrom=None,
                       ref_seq_len=0,
                       ref_seq_block_start=0,
                       ref_seq_block=None,
                       p=0.01,
                       start_snps_frac=0.0,
                       stop_snps_frac=1.0,
                       poisson_rng_seed=1,
                       base_sub_rng_seed=1,
                       prev_state=None,
                       **kwargs):
  """Given a chunk of the ref_seq get us some SNPs
  This is the stock generator and returns snp locations in a poisson distributed fashion.
  Inputs:
    chrom                - chromosome
    ref_seq_len          - length of whole sequence (needed to compute start and stop)
    ref_seq_block_start  - first coordinate of ref_seq_block
    ref_seq_block        - relevant chunk of the reference sequence
    p                    - probability of SNPs
    start_snps_frac      - start generating snps from here (0.0, 1.0)
    stop_snps_frac       - stop generating snps after this (0.0, 1.0) stop_snps_frac > start_snps_frac
    poisson_rng_seed     - SNP locator rng numpy.random.RandomState(seed)
    base_sub_rng_seed    - rng used to select ALT bases
    prev_state           - state of the model from the last call.
    kwargs               - absorbs any other parameters it does not use (and might pass on to sub functions)

  Outputs:
    variants       - (POS, REF, ALT)

  Notes: If the caller only calls this function when the ref_seq_block is at least partially within start_snps
  and stop_snps we will be faster, since otherwise this function simply returns (correctly) no variants and
  we needelessly use function overhead.

  Algorithm:
  1. Take the last SNP position
  2. Generate poisson distributed intervals
  3. Add them to the last SNP position
  4. Cut out any generated beyond the current block
  5. Using proper offset find the relevant REFs and generate ALTs, tossing out any unsuitable ones


  Test with one block
  >>> args = {'p': .1, 'poisson_rng_seed': 1, 'base_sub_rng_seed': 2}
  >>> ref_seq='ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT'
  >>>
  >>> variants, prev_state = candidate_variants(ref_seq_block_start=0, ref_seq_len=len(ref_seq), ref_seq_block=ref_seq, **args); print variants
  [(9, 'C', 'G'), (15, 'T', 'G'), (22, 'G', 'A'), (31, 'T', 'A'), (40, 'A', 'G'), (49, 'C', 'G'), (55, 'T', 'A')]

  Test with multiple blocks - should be same answer even though we have changed the blocks
  >>> variants = []
  >>> prev_state = None
  >>> for n,seq in enumerate(ref_seq):
  ...       _variants, prev_state = candidate_variants(ref_seq_block_start=n, ref_seq_len=len(ref_seq), ref_seq_block=seq, prev_state=prev_state, **args)
  ...       if len(_variants): variants += _variants
  ...
  >>> print variants
  [(9, 'C', 'G'), (15, 'T', 'G'), (22, 'G', 'A'), (31, 'T', 'A'), (40, 'A', 'G'), (49, 'C', 'G'), (55, 'T', 'A')]

  Test with offset.
  >>> args = {'p': 1, 'start_snps_frac': 0.33 , 'stop_snps_frac': 0.5, 'poisson_rng_seed': 1, 'base_sub_rng_seed': 2}
  >>> variants, prev_state = candidate_variants(ref_seq_block_start=0, ref_seq_len=len(ref_seq), ref_seq_block=ref_seq, **args); print variants
  [(20, 'A', 'T'), (21, 'C', 'A'), (22, 'G', 'A'), (23, 'T', 'A'), (24, 'A', 'G'), (25, 'C', 'G'), (26, 'G', 'C'), (27, 'T', 'G'), (28, 'A', 'C'), (29, 'C', 'T')]

  Test with offset and blocks, answer should be the same
  >>> variants = []
  >>> prev_state = None
  >>> for n,seq in enumerate(ref_seq):
  ...       _variants, prev_state = candidate_variants(ref_seq_block_start=n, ref_seq_len=len(ref_seq), ref_seq_block=seq, prev_state=prev_state, **args)
  ...       if len(_variants): variants += _variants
  ...
  >>> print variants
  [(20, 'A', 'T'), (21, 'C', 'A'), (22, 'G', 'A'), (23, 'T', 'A'), (24, 'A', 'G'), (25, 'C', 'G'), (26, 'G', 'C'), (27, 'T', 'G'), (28, 'A', 'C'), (29, 'C', 'T')]
  """
  start_snps = int(start_snps_frac * ref_seq_len)
  stop_snps = int(stop_snps_frac * ref_seq_len)
  try_to_run = True
  len_ref_seq = len(ref_seq_block)
  if ref_seq_block_start > stop_snps: try_to_run = False
  ref_seq_block_end = ref_seq_block_start + len(ref_seq_block)
  if ref_seq_block_end < start_snps: try_to_run = False
  poisson_blk_size = max(1, int(len_ref_seq * p * 1.01))  # Good guess as to how many SNPs in this interval
  stop_loc = min(ref_seq_block_end, stop_snps)

  if try_to_run:  # Only do this once we are ready to start generating SNPs
    keep_looping = True
    if prev_state is None:  # If we are running this for the first time, we need to initialize the base sub rng
      poisson_rng = numpy.random.RandomState(seed=poisson_rng_seed)
      base_sub_rng = numpy.random.RandomState(seed=base_sub_rng_seed)
      next_snp_locs = None
    else:
      poisson_rng = prev_state['poisson_rng']
      base_sub_rng = prev_state['base_sub_rng']
      next_snp_locs = prev_state['next_snp_locs']
  else:
    keep_looping = False

  variants = []
  while keep_looping:
    poiss = numpy.maximum(poisson_rng.poisson(lam=1.0 / p, size=poisson_blk_size),1)  # Need to move by 1 at least
    if next_snp_locs is None:  # These are the first SNPs we are generating
      these_snp_locs = numpy.cumsum(poiss) + start_snps
      # Our calculations start wrt the start of where we want to make SNPs. If we have gotten this far it means
      # that our ref_seq_block starts before start_snps and so we can do this computation
    else:
      these_snp_locs = numpy.concatenate((next_snp_locs, numpy.cumsum(poiss) + next_snp_locs[-1]))
      # We first try and use up the SNP locations we have already generated. This ensures that our simulations
      # are consistent even if we change block size
    if these_snp_locs[-1] > stop_loc:  # We generated more than we need.
      last_idx = numpy.argmax(these_snp_locs >= stop_loc)
      keep_looping = False
    else:
      last_idx = -1

    next_snp_locs = these_snp_locs[last_idx:]  # Will be used in the next run through, either this call or a future call
    these_snp_locs = these_snp_locs[:last_idx]

    base_subs = base_sub_rng.randint(3, size=len(these_snp_locs))
    # TODO: rewrite as list comprehension
    these_variants = []
    for loc, bs in zip(these_snp_locs - ref_seq_block_start, base_subs):
      # We need the relative positions since ref_seq is a substr
      if loc < 0: continue  # Before our window
      ref = ref_seq_block[loc]
      if ref not in ['G', 'A', 'T', 'C']: continue  # Not a valid base
      these_variants += [(loc + ref_seq_block_start, ref, base_sub_mat[ref][bs])]

    if len(these_variants) > 0:
      variants += these_variants

  if try_to_run:  # This flag is used to avoid exiting the function at multiple points.
    prev_state = {
      'poisson_rng': poisson_rng,
      'base_sub_rng': base_sub_rng,
      'next_snp_locs': next_snp_locs
    }

  return variants, prev_state


if __name__ == "__main__":
  import doctest
  doctest.testmod()

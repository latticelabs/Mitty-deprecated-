"""TODO: Write docs"""
import numpy
import logging

logger = logging.getLogger(__name__)


def candidate_variants(chrom=None, start_loc=0, ref_seq=None, prev_state=None,
                       p=0.01, poisson_rng_seed=1, base_sub_rng_seed=1):
  """Given a chunk of the ref_seq get us some SNPs
  This is the stock generator and returns snp locations in a poisson distributed fashion.
  Inputs:
    chrom          - chromosome
    start_loc      - first coordinate of ref_seq
    ref_seq        - relevant chunk of the reference sequence
    rng            - numpy.random.RandomState(seed)
    prev_state     - state of the model from the last call.
                     We store the base substit rng here too
                     Leave at default for first run.
                     Stores location of the next SNP otherwise
    p              - probability of SNP at any given location (passed as part of kwargs)

  Outputs:
    variants       - (POS, REF, ALT)

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
  >>>
  >>> variants, prev_state = candidate_variants(start_loc=0, ref_seq=ref_seq, **args)
  >>> print variants
  [(9, 'C', 'G'), (15, 'T', 'G'), (22, 'G', 'A'), (31, 'T', 'A'), (40, 'A', 'G'), (49, 'C', 'G'), (55, 'T', 'A')]

  Test with multiple blocks - should be same answer even though we have changed the blocks
  >>> variants = []
  >>> prev_state = None
  >>> for n,seq in enumerate(ref_seq):
  ...       _variants, prev_state = candidate_variants(start_loc=n, ref_seq=seq, prev_state=prev_state, **args)
  ...       if len(_variants): variants += _variants
  ...
  >>> print variants
  [(9, 'C', 'G'), (15, 'T', 'G'), (22, 'G', 'A'), (31, 'T', 'A'), (40, 'A', 'G'), (49, 'C', 'G'), (55, 'T', 'A')]
  """
  base_sub_mat = {  # GATC
                    'G': 'ATC',
                    'A': 'TCG',
                    'T': 'CGA',
                    'C': 'GAT'
  }

  variants = []
  len_ref_seq = len(ref_seq)
  stop_loc = start_loc + len_ref_seq
  if prev_state is None:  # If we are running this for the first time, we need to initialize the base sub rng
    poisson_rng = numpy.random.RandomState(seed=poisson_rng_seed)
    base_sub_rng = numpy.random.RandomState(seed=base_sub_rng_seed)
    next_snp_locs = None
  else:
    poisson_rng = prev_state['poisson_rng']
    base_sub_rng = prev_state['base_sub_rng']
    next_snp_locs = prev_state['next_snp_locs']

  poisson_blk_size = max(1, int(len_ref_seq * p * 1.01))  # Good guess as to how many SNPs in this interval

  keep_looping = True
  while keep_looping:
    poiss = poisson_rng.poisson(lam=1.0 / p, size=poisson_blk_size)
    if next_snp_locs is None:
      these_snp_locs = numpy.cumsum(poiss) + start_loc  # Our calculations start wrt the start of the block we are given
    else:
      these_snp_locs = numpy.concatenate((next_snp_locs, numpy.cumsum(poiss) + next_snp_locs[-1]))
      # Our calculations start wrt the next SNP with the next SNP being the first in line

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
    for loc, bs in zip(these_snp_locs - start_loc, base_subs):
      # We need the relative positions since ref_seq is a substr
      if loc < 0: continue  # Before our window
      ref = ref_seq[loc]
      if ref not in ['G', 'A', 'T', 'C']: continue  # Not a valid base
      these_variants += [(loc + start_loc, ref, base_sub_mat[ref][bs])]

    if len(these_variants) > 0:
      variants += these_variants

  prev_state = {
    'poisson_rng': poisson_rng,
    'base_sub_rng': base_sub_rng,
    'next_snp_locs': next_snp_locs
  }

  return variants, prev_state


if __name__ == "__main__":
  import doctest

  doctest.testmod()

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
het_type = ['0/1', '1/0']  # The two types of het


def variant(ref_seq=None, ref_seq_len=0,
            phet=0.0, p=0.01,
            start_snps_frac=0.0, stop_snps_frac=1.0,
            het_rng_seed=1,
            strand_rng_seed=4,
            poisson_rng_seed=2,
            base_sub_rng_seed=3,
            block_size=10000, **kwargs):
  """A generator which returns a variant when asked for. This is the stock SNP generator and returns snp locations in
  a poisson distributed fashion.
  Inputs:
    ref_seq              - The reference sequence
    ref_seq_len          - length of whole sequence (needed to compute start and stop)
    phet                 - probability of having heterozygous reads
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
  def get_locs_and_subs():  # Simply a convenience.
    het_or_not = het_rng.rand(block_size)
    strand_no = strand_rng.randint(2, size=block_size)
    poiss = numpy.maximum(poisson_rng.poisson(lam=1.0 / p, size=block_size), 1)  # Please sir, can I have some more?
    base_subs = base_sub_rng.randint(3, size=block_size)
    locs = numpy.cumsum(poiss) + start_offset
    return het_or_not, strand_no, locs, base_subs, locs[-1]

  het_rng = numpy.random.RandomState(seed=het_rng_seed)
  strand_rng = numpy.random.RandomState(seed=strand_rng_seed)
  poisson_rng = numpy.random.RandomState(seed=poisson_rng_seed)
  base_sub_rng = numpy.random.RandomState(seed=base_sub_rng_seed)

  start_offset = int(ref_seq_len * start_snps_frac)
  snp_end = int(ref_seq_len * stop_snps_frac)

  het_or_not, strand_no, locs, subs, start_offset = get_locs_and_subs()
  internal_cntr = 0
  while locs[internal_cntr] < snp_end:
    vl = locs[internal_cntr]
    ref = ref_seq[vl]
    if ref in ['G', 'A', 'T', 'C']:
      alt = base_sub_mat[ref][subs[internal_cntr]]
      gt = '1/1' if het_or_not[internal_cntr] > phet else het_type[strand_no[internal_cntr]]
    else:
      alt = None  # Not valid, skip to next

    internal_cntr += 1
    if internal_cntr == locs.size:  # Ran out of random numbers, need to generate some more
      het_or_not, strand_no, locs, subs, start_offset = get_locs_and_subs()
      internal_cntr = 0

    if alt is not None:
      yield (vl, ref, alt, gt, vl + 2, None)  # POS, REF, ALT, skipto, list(footprints)
                                          # footprints, in this case, is None, since we simply skip forward
      # We have vl + 2 because we want 1 base buffer between variants, even SNPs (see Readme)

if __name__ == "__main__":
  import doctest
  doctest.testmod()

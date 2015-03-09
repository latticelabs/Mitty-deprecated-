"""This is the stock SNP plugin. It uses four independent RNGs to locate SNPs along a reference sequence, assign each
SNP a zygosity and assign an ALT base.
"""
import mitty.lib.util as mutil
import logging
logger = logging.getLogger(__name__)

__example_param_text = """
{
  "chromosome": [4],      # List of chromosomes to apply the data to
  "p": 0.01,              # probability that the SNP will happen at any given base
  "phet": 0.5,            # probability that the data will be heterozygous
  "t_mat": [[ 0.32654629,  0.17292732,  0.24524503,  0.25528135],  # Base transition matrix
            [ 0.3489394,   0.25942695,  0.04942584,  0.3422078],   # Leading diagonal is ignored
            [ 0.28778188,  0.21087004,  0.25963262,  0.24171546],
            [ 0.21644706,  0.20588717,  0.24978216,  0.32788362]],
}
"""

_description = """
This is the stock SNP plugin. A typical parameter set resembles
""" + __example_param_text

_example_params = eval(__example_param_text)


def variant_generator(ref={},
                      chromosome=None,
                      phet=0.5,
                      p=0.01,
                      t_mat=None,
                      master_seed=1,
                      **kwargs):
  assert 0 <= p <= 1.0, "Probability out of range"
  assert 0 <= phet <= 1.0, "Probability out of range"
  if t_mat is None:
    t_mat = [[0.32654629, 0.17292732, 0.24524503, 0.25528135],
             [0.3489394, 0.25942695, 0.04942584, 0.3422078],
             [0.28778188, 0.21087004, 0.25963262, 0.24171546],
             [0.21644706, 0.20588717, 0.24978216, 0.32788362]]

  logger.debug('Master seed: {:d}'.format(master_seed))
  base_loc_rng, base_t_rng, het_rng, copy_rng = mutil.initialize_rngs(master_seed, 4)

  for chrom in chromosome:
    ref_chrom = ref[chrom]
    snp_locs = [x for x in mutil.place_poisson(base_loc_rng, p, len(ref_chrom)) if ref_chrom[x] != 'N']
    het_type = mutil.zygosity(len(snp_locs), phet, het_rng, copy_rng)
    base_subs = mutil.base_subs(ref_chrom, snp_locs, t_mat, base_t_rng)
    yield {chrom: [[pos + 1 for pos in snp_locs], [ref_chrom[pos] for pos in snp_locs], base_subs, het_type]}
    # +1 because VCF files are 1 indexed
    #alts will be 0 if ref is not one of ACTG


def test():
  """Basic test"""
  ref = {
    1: 'ACTGACTGACTG',
    2: 'ACTGACTGACTGACTGACTGACTG'
  }
  vg = variant_generator(ref, chromosome=[1, 2], p=1.0)
  for v_list in vg:
    for chrom, dnv in v_list.iteritems():
      for p, r, a, gt in zip(*dnv):
        assert ref[chrom][p - 1] == r


if __name__ == "__main__":
  print _description
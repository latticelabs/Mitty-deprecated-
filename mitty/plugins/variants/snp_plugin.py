"""This is the stock SNP plugin. It uses four independent RNGs to locate SNPs along a reference sequence, assign each
SNP a zygosity and assign an ALT base.
"""
import numpy as np

import mitty.lib
import mitty.lib.util as mutil
import logging
logger = logging.getLogger(__name__)

__example_param_text = """
{
  "p": 0.01,              # probability that the SNP will happen at any given base
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


class Model:
  def __init__(self, ref=[], p=0.01, t_mat=None):
    assert 0 <= p <= 1.0, "Probability out of range"
    if t_mat is None:
      t_mat = [[0.32654629, 0.17292732, 0.24524503, 0.25528135],
               [0.3489394, 0.25942695, 0.04942584, 0.3422078],
               [0.28778188, 0.21087004, 0.25963262, 0.24171546],
               [0.21644706, 0.20588717, 0.24978216, 0.32788362]]
    self.p, self.t_mat = p, t_mat

  def get_variants(self, ref, chrom, p, f, seed=1):
    """This function is called by the simulator to obtain variants.

    :param ref: reference sequence as a string
    :param chrom: chromosome number (1,2,3,4...)
    :param p: array/list of probability values
    :param f: array/list of frequency values
    :param seed: seed for the random number generators
    :return: 5 arrays/lists/iterables all of the same length
              pos   - position of SNPs
              stop  - stop locations, (pos + 1 for SNPs)
              ref   - reference base,
              alt   - alt base,
              p     - probability value for this variant. These are uniformly distributed random values
    """
    assert 0 < seed < mitty.lib.SEED_MAX
    logger.debug('Master seed: {:d}'.format(seed))
    assert type(p) == np.ndarray, 'Site Freq spectrum should be numpy array'
    assert type(f) == np.ndarray, 'Site Freq spectrum should be numpy array'
    assert abs(1.0 - f.sum()) < 1e-6, 'Site freq spectrum should sum to 1.0'

    p_eff = self.p / (p * f).sum()  # Effective per-base probability for master list

    base_loc_rng, base_t_rng, freq_rng = mutil.initialize_rngs(seed, 3)
    snp_locs = np.array([x for x in mutil.place_poisson(base_loc_rng, p_eff, len(ref)) if ref[x] != 'N'], dtype='i4')
    base_subs = mutil.base_subs(ref, snp_locs, self.t_mat, base_t_rng)

    # +1 because VCF files are 1 indexed
    #alts will be 0 if ref is not one of ACTG
    return snp_locs + 1, snp_locs + 2, [ref[n] for n in snp_locs], base_subs, freq_rng.rand(len(snp_locs))


def test():
  """Basic test"""
  ref_seq = 'ACTGACTGACTGACTGACTGACTGACTGACTGACTG'
  m = Model(p=0.1)
  pos, stop, ref, alt, p = m.get_variants(ref_seq[1], 1, np.array([0.2]), np.array([1.0]), seed=10)
  for p, r in zip(pos, ref):
    assert r == ref_seq[pos - 1]


if __name__ == "__main__":
  print _description
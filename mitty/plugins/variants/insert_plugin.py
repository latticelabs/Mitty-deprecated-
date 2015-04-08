"""This is the stock insertion generator"""
import numpy as np

import mitty.lib
import mitty.lib.util as mutil
from mitty.plugins.variants import scale_probability_and_validate

import logging
logger = logging.getLogger(__name__)

__example_param_text = """
{
  "p": 0.0001,        # Per-base probability of having an insertion
  "t_mat": [[ 0.32654629,  0.17292732,  0.24524503,  0.25528135],  # Base transition matrix
            [ 0.3489394,   0.25942695,  0.04942584,  0.3422078],
            [ 0.28778188,  0.21087004,  0.25963262,  0.24171546],
            [ 0.21644706,  0.20588717,  0.24978216,  0.32788362]],
  "p_end": 0.1,       # Probability of chain ending
  "max_len": 1000     # Maximum length of insertion
}
"""

_description = """
Stock insertion model that generates sequences with same base transition matrix as the human genome and creates a
power-law distribution of insertion lengths.
A typical parameter set resembles
""" + __example_param_text

_example_params = eval(__example_param_text)


class Model:
  def __init__(self, p=0.01, t_mat=None, p_end=0.1, max_len=1000, **kwargs):
    assert 0 <= p <= 1.0, "Probability out of range"
    assert 0 <= p_end <= 1.0, "Probability out of range"
    assert 0 < max_len, 'max_len needs to be 1 or more'

    if t_mat is None:
      t_mat = [[0.32654629, 0.17292732, 0.24524503, 0.25528135],
               [0.3489394, 0.25942695, 0.04942584, 0.3422078],
               [0.28778188, 0.21087004, 0.25963262, 0.24171546],
               [0.21644706, 0.20588717, 0.24978216, 0.32788362]]
    self.p, self.t_mat, self.p_end, self.max_len = p, t_mat, p_end, max_len

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

    base_loc_rng, ins_markov_rng = mutil.initialize_rngs(seed, 2)

    pt_mat = mutil.add_p_end_to_t_mat(self.t_mat, self.p_end)
    p_eff = scale_probability_and_validate(self.p, p, f)
    ins_locs = mutil.place_poisson_seq(base_loc_rng, p_eff, 0, len(ref), ref)  #np.array([x for x in mutil.place_poisson(base_loc_rng, p_eff, 0, len(ref)) if ref[x] != 'N'], dtype='i4')
    ins_list, len_list = mutil.markov_sequences(ref, ins_locs, self.max_len, pt_mat, ins_markov_rng)
    lengths = np.array(len_list, dtype='i4')

    return ins_locs, ins_locs + 1, [ins[0] for ins in ins_list], ins_list, lengths / float(lengths.max())


def test():
  """Basic test"""
  ref_seq = 'ACTGACTGACTGACTGACTGACTGACTGACTGACTG'
  m = Model(p=0.1)
  pos, stop, ref, alt, p = m.get_variants(ref_seq, 1, np.array([0.2]), np.array([1.0]), seed=10)
  for p, r in zip(pos, alt):
    assert r[0] == ref_seq[p]


if __name__ == "__main__":
    print _description
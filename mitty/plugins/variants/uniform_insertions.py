"""This insertion generator generates insertions with uniformly distributed lengths"""
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
  "min_len": 2,    # Minimum length of insertion
  "max_len": 30    # Maximum length of insertion
}
"""

_description = """
Insertion model that generates sequences with same base transition matrix as the human genome and creates a
uniform distribution of insertion lengths. A typical parameter set resembles
""" + __example_param_text

_example_params = eval(__example_param_text)


class Model:
  def __init__(self, p=0.01, t_mat=None, min_len=2, max_len=30, **kwargs):
    assert 0 <= p <= 1.0, "Probability out of range"
    assert 0 < min_len, 'min_len needs to be 1 or more'
    assert min_len <= max_len, 'max_len needs to be >= than min_len'

    if t_mat is None:
      t_mat = [[0.32654629, 0.17292732, 0.24524503, 0.25528135],
               [0.3489394, 0.25942695, 0.04942584, 0.3422078],
               [0.28778188, 0.21087004, 0.25963262, 0.24171546],
               [0.21644706, 0.20588717, 0.24978216, 0.32788362]]
    self.p, self.t_mat, self.min_len, self.max_len = p, t_mat, min_len, max_len

  def get_variants(self, ref, p=None, f=None, seed=1, **kwargs):
    """This function is called by the simulator to obtain variants.

    :param ref: reference sequence as a string
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

    base_loc_rng, ins_markov_rng, ins_len_rng = mutil.initialize_rngs(seed, 3)

    pt_mat = mutil.add_p_end_to_t_mat(self.t_mat, 0.0)  # The sequence will only end at max len
    p_eff = scale_probability_and_validate(self.p, p, f)
    ins_locs = mutil.place_poisson_seq(base_loc_rng, p_eff, 0, len(ref), ref)  #np.array([x for x in mutil.place_poisson(base_loc_rng, p_eff, 0, len(ref)) if ref[x] != 'N'], dtype='i4')
    ins_lens = ins_len_rng.randint(low=self.min_len, high=self.max_len + 1, size=len(ins_locs))
    ins_list, len_list = mutil.markov_sequences(ref, ins_locs, ins_lens, pt_mat, ins_markov_rng)
    lengths = np.array(len_list, dtype='i4')

    return ins_locs, ins_locs + 1, [ins[0] for ins in ins_list], ins_list, (1.0 - lengths / float(lengths.max())) if lengths.shape[0] else []


def test0():
  """Edge case - no variants generated"""
  ref_seq = 'ACTGACTGACTGACTGACTGACTGACTGACTGACTG'
  m = Model(p=0.00001)
  pos, stop, ref, alt, p = m.get_variants(ref_seq, seed=10)
  assert len(pos) == 0  # This should just run and not crash


def test():
  """Basic test"""
  ref_seq = 'ACTGACTGACTGACTGACTGACTGACTGACTGACTG'
  m = Model(p=0.1)
  pos, stop, ref, alt, p = m.get_variants(ref_seq, seed=10)
  for p, r in zip(pos, alt):
    assert r[0] == ref_seq[p]


if __name__ == "__main__":
    print _description
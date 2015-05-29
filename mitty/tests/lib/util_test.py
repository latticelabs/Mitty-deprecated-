import numpy
from numpy.testing import assert_array_equal, assert_array_almost_equal

import mitty.lib.util


def place_poisson_test():
  """Poisson RNG"""
  rng = numpy.random.RandomState(seed=1)
  p = .01
  end_p = 10000
  correct_locs = rng.geometric(p=p, size=end_p).cumsum()
  idx, = numpy.nonzero(correct_locs >= end_p)
  correct_locs = correct_locs[:idx[0]]

  rng = numpy.random.RandomState(seed=1)
  computed_locs = mitty.lib.util.place_poisson_seq(rng, p, 0, end_p, 'A' * end_p)

  assert_array_equal(correct_locs, computed_locs)


class MockRng:
  def __init__(self, r):
    self.r = r

  def rand(self, n):
    return numpy.array(self.r[:n], dtype=numpy.double)


def base_subs_test():
  """Base substitution."""

  seq = 'ACGT'
  sub_pts = [0, 1, 2, 3]
  t_mat = [[0.0, 0.3, 0.3, 0.3],
           [0.3, 0.0, 0.3, 0.3],
           [0.3, 0.3, 0.0, 0.3],
           [0.3, 0.3, 0.3, 0.0]]
  rng = MockRng([0.5, 0.5, 0.5, 0.5])
  # -> G G C C

  base_subbed = mitty.lib.util.base_subs(seq, sub_pts, t_mat, rng)
  assert ['G', 'G', 'C', 'C'] == base_subbed, base_subbed


def add_p_end_to_t_mat_test():
  """Adding p_end to t_mat."""
  t_mat = [[0.25, 0.25, 0.25, 0.25],
           [0.1, 0.1, 0.1, 0.1],
           [0.2, 0.1, 0.2, 0.1],
           [0.1, 0.2, 0.2, 0.1]]
  p_end = 0.1
  correct_t_mat = [[0.225, 0.225, 0.225, 0.225, 0.1],
                   [0.225, 0.225, 0.225, 0.225, 0.1],
                   [0.3, 0.15, 0.3, 0.15, 0.1],
                   [0.15, 0.3, 0.3, 0.15, 0.1]]
  assert_array_almost_equal(correct_t_mat, mitty.lib.util.add_p_end_to_t_mat(t_mat, p_end))


def sequence_gen_test():
  """Markov chain sequence generator, normal termination."""
  seq = 'A'
  ins_pts = [0]
  max_len = 10
  #          A    C    G    T    x
  t_mat = [[0.2, 0.2, 0.2, 0.2, 0.2],
           [0.1, 0.1, 0.1, 0.1, 0.6],
           [0.2, 0.1, 0.2, 0.1, 0.4],
           [0.1, 0.2, 0.2, 0.1, 0.4]]
  rng = MockRng([0.5, .72, .0001, .3, .147, .092, .186, .345, .397, 1.0])
  #   0.5, .72,        .0001,  .3,  .147, .092, .186, .345, .397, 1.0
  # A->G -> x(ignored)-> A->   C->   C->    A->   A->   C->   T->  (end)
  seq_l, l = mitty.lib.util.markov_sequences(seq, ins_pts, max_len, t_mat, rng)
  assert seq_l[0] == 'AGACCAACT', seq_l[0]
  assert l[0] == 9


def sequence_gen_test2():
  """Markov chain sequence generator, terminate when too long."""
  seq = 'A'
  ins_pts = [0]
  max_len = 7
  #          A    C    G    T    x
  t_mat = [[0.2, 0.2, 0.2, 0.2, 0.2],
           [0.1, 0.1, 0.1, 0.1, 0.6],
           [0.2, 0.1, 0.2, 0.1, 0.4],
           [0.1, 0.2, 0.2, 0.1, 0.4]]
  rng = MockRng([0.5, .72, .0001, .3, .147, .092, .186, .345, .397, 0.55])
  #   0.5, .72,        .0001,  .3,  .147, .092, .186, .345, .397,  0.55
  # A->G -> x(ignored)-> A->   C->   C->    A->   A->  (force end)
  seq_l, l = mitty.lib.util.markov_sequences(seq, ins_pts, max_len, t_mat, rng)
  assert seq_l[0] == 'AGACCAA', seq_l[0]
  assert l[0] == 7
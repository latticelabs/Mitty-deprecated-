import numpy
from numpy.testing import assert_array_equal, assert_array_almost_equal

import mitty.lib.util


def place_poisson_test():
  rng = numpy.random.RandomState(seed=1)
  p = .01
  end_p = 10000
  correct_locs = rng.poisson(lam=1./p, size=end_p).cumsum()
  idx, = numpy.nonzero(correct_locs >= end_p)
  correct_locs = correct_locs[:idx[0]]

  rng = numpy.random.RandomState(seed=1)
  computed_locs = mitty.lib.util.place_poisson(rng, p, end_p)

  assert_array_equal(correct_locs, computed_locs)


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
  class MockRng:
    def __init__(self, r):
      self.rand_iter = (n for n in r)

    def rand(self):
      return self.rand_iter.next()

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
  seq_l = mitty.lib.util.markov_sequences(seq, ins_pts, max_len, t_mat, rng)
  assert seq_l[0] == 'GACCAACT', seq_l[0]


def sequence_gen_test2():
  """Markov chain sequence generator, terminate when too long."""
  class MockRng:
    def __init__(self, r):
      self.rand_iter = (n for n in r)

    def rand(self):
      return self.rand_iter.next()

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
  # A->G -> x(ignored)-> A->   C->   C->    A->   A->   C-> (force end)
  seq_l = mitty.lib.util.markov_sequences(seq, ins_pts, max_len, t_mat, rng)
  assert seq_l[0] == 'GACCAAC', seq_l[0]
import mitty.lib.variants as vr
from nose.tools import assert_sequence_equal
from numpy.testing import assert_array_equal


def initialize_test():
  """Initialize variant list"""
  pos = [1, 10, 20]
  stop = [2, 11, 21]
  ref = ['A', 'C', 'T']
  alt = ['AA', 'CAT', 'G']
  p = [0.1, 0.5, 0.9]
  l = vr.VariantList(pos, stop, ref, alt, p)
  assert l.variants['pos'][1] == 10
  assert l.variants['alt'][2] == 'G'


def add_test():
  """Add variations to the list."""
  pos = [1, 10, 20]
  stop = [2, 11, 21]
  ref = ['A', 'C', 'T']
  alt = ['AA', 'CAT', 'G']
  p = [0.1, 0.5, 0.9]
  l = vr.VariantList(pos, stop, ref, alt, p)

  pos = [5, 15]
  stop = [6, 16]
  ref = ['T', 'G']
  alt = ['A', 'GAT']
  p = [0.1, 0.5]
  l.add(pos, stop, ref, alt, p)

  assert l.variants['pos'][1] == 10
  assert l.variants['alt'][2] == 'G'

  assert l.variants['pos'][4] == 15
  assert l.variants['alt'][3] == 'A'


def sort_test():
  """Sort variations"""
  pos = [1, 10, 20]
  stop = [2, 11, 21]
  ref = ['A', 'C', 'T']
  alt = ['AA', 'CAT', 'G']
  p = [0.1, 0.5, 0.9]
  l = vr.VariantList(pos, stop, ref, alt, p)

  pos = [5, 15]
  stop = [6, 16]
  ref = ['T', 'G']
  alt = ['A', 'GAT']
  p = [0.1, 0.5]
  l.add(pos, stop, ref, alt, p)

  l.sort()

  assert l.variants['pos'][1] == 5
  assert l.variants['alt'][2] == 'CAT'


def select_test():
  """Select variants based on probability"""
  pos = [1, 10, 20]
  stop = [2, 11, 21]
  ref = ['A', 'C', 'T']
  alt = ['AA', 'CAT', 'G']
  p = [0.1, 0.5, 0.9]
  l = vr.VariantList(pos, stop, ref, alt, p)

  class MyRng:
    def __init__(self):
      self.v = vr.np.array([[0.1, 0.4, 0.4], [0.05, 0.8, 0.1]]).T

    def rand(self, d0, d1):
      return self.v

  rng = MyRng()
  idx = l.select(rng)
  assert_array_equal(idx[0], [1, 2])
  assert_array_equal(idx[1], [0, 2])


def zip_test():
  """Zip chromosomes together (Need to make tests more rigorous)."""
  pos = [1, 2, 20]
  stop = [5, 3, 21]
  ref = ['ACTGA', 'C', 'T']
  alt = ['A', 'CAT', 'G']
  p = [0.1, 0.5, 0.9]
  l = vr.VariantList(pos, stop, ref, alt, p)
  chrom = l.zip_up_chromosome([0, 1, 2], [0, 1, 2])
  assert_sequence_equal(chrom, [(0, 2), (2, 2)])



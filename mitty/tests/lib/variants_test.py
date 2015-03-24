import mitty.lib.variants as vr
from nose.tools import assert_sequence_equal
from numpy.testing import assert_array_equal, assert_array_almost_equal


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


def balance_test():
  """Balance the site frequency spectrum"""
  pos = [1, 10, 20, 30]
  stop = [2, 11, 21, 31]
  ref = ['A', 'C', 'T', 'G']
  alt = ['T', 'G', 'G', 'A']
  p = [0.1, 0.9, 0.3, 0.2]
  l = vr.VariantList(pos, stop, ref, alt, p)

  # The core things here are to make sure the sfs is balanced AND the rank ordering is maintained as far as possible
  p_sfs = [0.1, 0.2, 0.3]
  f_sfs = [0.5, 0.25, 0.25]

  l.balance_probabilities(p_sfs, f_sfs)
  assert_array_almost_equal(l.variants['p'], [0.1, 0.3, 0.2, 0.1], decimal=3)


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


def collision_test():
  """Collisions"""
  pos =  [1, 3, 5, 7]
  stop = [3, 4, 7, 8]
  #          ^     ^  <-- these collide
  z_idx = vr.avoid_collisions(pos, stop, [0, 1, 2, 3])
  assert_array_equal(z_idx, [0, 2])


def merge_test():
  """Merging"""
  pos =  [1, 3, 5, 7]
  z0 = [0, 2, 3]
  z1 = [1, 2, 3]
  chrom = vr.merge_homozygous(pos, z0, z1)
  assert_array_equal(chrom, [(0, 0), (1, 1), (2, 2), (3, 2)])


def zip_test():
  """Zip chromosomes together"""
  pos = [1, 2, 20]
  stop = [5, 3, 21]
  ref = ['ACTGA', 'C', 'T']
  alt = ['A', 'CAT', 'G']
  p = [0.1, 0.5, 0.9]
  l = vr.VariantList(pos, stop, ref, alt, p)
  chrom = l.zip_up_chromosome([0, 1, 2], [0, 1, 2])
  assert_sequence_equal(chrom, [(0, 2), (2, 2)])



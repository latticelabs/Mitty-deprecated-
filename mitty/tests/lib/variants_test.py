import mitty.lib.variants as vr

from nose.tools import assert_sequence_equal
from nose.tools import assert_raises

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal


def assert_index_array_equal(x, y):
  dt = [('index', 'i4'), ('gt', 'i1')]
  assert_array_equal(np.array(x, dtype=dt), np.array(y, dtype=dt))


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
  assert_index_array_equal(chrom, [(0, 0), (1, 1), (2, 2), (3, 2)])


def merge_test1():
  """Merge stragglers"""
  pos = [1, 2, 3, 4, 5]
  z0 = [2, 4]
  z1 = [3]
  chrom = vr.merge_homozygous(pos, z0, z1)
  assert_index_array_equal(chrom, [(2, 0), (3, 1), (4, 0)])

  pos = [1, 2, 3, 4, 5]
  z0 = [0, 1]
  z1 = [2, 3]
  chrom = vr.merge_homozygous(pos, z0, z1)
  assert_index_array_equal(chrom, [(0, 0), (1, 0), (2, 1), (3, 1)])


def merge_test2():
  """Merging staggered"""
  pos = [1, 2, 3, 4, 5]
  z0 = [1, 2, 4]
  z1 = [0, 3]
  chrom = vr.merge_homozygous(pos, z0, z1)
  assert_index_array_equal(chrom, [(0, 1), (1, 0), (2, 0), (3, 1), (4, 0)])


def zip_test():
  """Zip chromosomes together"""
  pos = [1, 2, 20]
  stop = [5, 3, 21]
  ref = ['ACTGA', 'C', 'T']
  alt = ['A', 'CAT', 'G']
  p = [0.1, 0.5, 0.9]
  l = vr.VariantList(pos, stop, ref, alt, p)
  chrom = l.zip_up_chromosome([0, 1, 2], [0, 1, 2])
  assert_index_array_equal(chrom, [(0, 2), (2, 2)])


def master_list_roundtrip_test():
  """Master list round trip"""
  genome_metadata = [{'seq_id': 'chr1', 'seq_len': 10, 'seq_md5': '10'}]

  pos = [1, 10, 20]
  stop = [2, 11, 21]
  ref = ['A', 'C', 'T']
  alt = ['AA', 'CAT', 'G']
  p = [0.1, 0.5, 0.9]
  ml = vr.VariantList(pos, stop, ref, alt, p)

  assert_raises(RuntimeError, vr.Population)  # We should complain about not having metadata since this is a new file

  pl = vr.Population(genome_metadata=genome_metadata)
  assert_raises(AssertionError, pl.set_master_list, 1, ml)  # We gotta sort this first

  ml.sort()
  pl.set_master_list(1, ml)

  ml2 = pl.get_master_list(1)
  assert_array_equal(ml.variants, ml2.variants)

  pl = vr.Population(genome_metadata=genome_metadata)
  pl.set_master_list(1, ml)
  assert_array_equal(ml.variants, pl.get_master_list(1).variants)

  from mitty.version import __version__
  assert pl.get_version() == __version__


def sample_roundtrip_test():
  """Sample round-trip (save and load)"""
  genome_metadata = [{'seq_id': 'chr1', 'seq_len': 10, 'seq_md5': '10'}]
  chrom = np.array([(1, 0), (2, 1), (3, 2), (1073741823, 2)], dtype=[('index', 'i4'), ('gt', 'i1')])
  pl = vr.Population(genome_metadata=genome_metadata)
  pl.add_sample_chromosome(chrom=1, sample_name='my_sample', indexes=chrom)
  c2 = pl.get_sample_chromosome(chrom=1, sample_name='my_sample')
  assert_array_equal(chrom, c2)


def chrom_metadata_roundtrip_test():
  """Chromosome metadata round-trip"""
  genome_metadata = [
    {'seq_id': 'Old McDonald had a farm', 'seq_len': 100, 'seq_md5': '23'},
    {'seq_id': 'Five little monkeys jumping on the bed', 'seq_len': 200, 'seq_md5': '99'},
  ]

  pl = vr.Population(genome_metadata=genome_metadata)
  assert_sequence_equal(genome_metadata, pl.get_genome_metadata())
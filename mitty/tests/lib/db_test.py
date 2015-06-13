from nose.tools import assert_sequence_equal
from nose.tools import assert_raises
from numpy.testing import assert_array_equal

import mitty.lib.variants as vr
import mitty.lib.db as mdb


def master_list_roundtrip_test():
  """Master list round trip"""
  pos = [1, 10, 20]
  stop = [2, 11, 21]
  ref = ['A', 'C', 'T']
  alt = ['AA', 'CAT', 'G']
  p = [0.1, 0.5, 0.9]
  ml = vr.VariantList(pos, stop, ref, alt, p)

  conn = mdb.connect(':memory:')
  assert_raises(AssertionError, mdb.save_master_list, conn, 1, ml)  # We gotta sort this first

  ml.sort()
  mdb.save_master_list(conn, 1, ml)

  ml2 = mdb.load_master_list(conn, 1)
  assert_array_equal(ml.variants, ml2.variants)


def sample_roundtrip_test():
  """Sample round-trip (save and load from database)"""
  chrom = [(1, 0), (2, 1), (3, 2), (1073741823, 2)]
  conn = mdb.connect(':memory:')

  cid = mdb.save_sample(conn, 0, 0, 1, chrom)
  assert cid == 1

  c2 = mdb.load_sample(conn, 0, 0, 1)
  assert_sequence_equal(chrom, c2)


def chrom_metadata_roundtrip_test():
  """Chromosome metadata round-trip"""
  conn = mdb.connect(':memory:')

  seq_id1 = 'Old McDonald had a farm'
  seq1 = 'Eeya Eeya O. And on this farm he had a goat, Eeya Eeya O'
  mdb.save_chromosome_metadata(conn, 1, seq_id1, len(seq1), '22')

  seq_id2 = 'Five little monkeys jumping on the bed'
  seq2 = 'One fell down and bumped his head. Mommy called the doctor, and the doctor said "No more monkeys jumping on the bed"'
  mdb.save_chromosome_metadata(conn, 2, seq_id2, len(seq2), '9675')

  _, sid1, slen1, md51 = mdb.load_chromosome_metadata(conn, 1)
  _, sid2, slen2, md52 = mdb.load_chromosome_metadata(conn, 2)

  assert sid1 == seq_id1
  assert sid2 == seq_id2
  assert len(seq1) == slen1
  assert len(seq2) == slen2
  assert md51 == '22'
  assert md52 == '9675'
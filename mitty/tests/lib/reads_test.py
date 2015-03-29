import mitty.lib.variants as vr
import mitty.lib.reads as reads
from nose.tools import assert_sequence_equal


def expand_seq_test1():
  """Sequence expand no variants"""
  ref_seq = 'ACTGACTGACTGACT'
  ml = vr.VariantList()
  chrom = []
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)

  assert ref_seq == alt_seq
  assert_sequence_equal(beacons, [])


def expand_seq_test2():
  """Sequence expand one variant of each type"""
  #          123456789012345
  ref_seq = 'ACTGACTGACTGACT'

  pos = [4]
  stop = [5]
  ref = ['G']
  alt = ['T']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 1)]
  #        123456789012345
  #ref     ACTGACTGACTGACT
  m_alt = 'ACTTACTGACTGACT'
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert ref_seq == alt_seq, alt_seq
  assert_sequence_equal(beacons, [], str(beacons))

  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 1)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons.tolist(), [(3, 3, 0)])

  pos = [4]
  stop = [5]
  ref = ['G']
  alt = ['GAAA']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 0)]

  #        1234   56789012345
  #ref     ACTG   ACTGACTGACT
  m_alt = 'ACTGAAAACTGACTGACT'
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons.tolist(), [(3, 3, 3)])

  pos = [4]
  stop = [8]
  ref = ['GACT']
  alt = ['G']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2)]

  #        123456789012345
  #            xxx
  #ref     ACTGACTGACTGACT
  m_alt = 'ACTGGACTGACT'
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons.tolist(), [(3, 3, -3)])


def expand_seq_test3():
  """Sequence expand multiple variants, different combinations"""
  #          123456789012
  ref_seq = 'ACTGACTGACTG'

  pos = [4, 6]
  stop = [5, 9]
  ref = ['G', 'CTG']
  alt = ['GAT', 'C']
  p = [0.1, 0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 1), (1, 1)]
  #        1234  5678   9012
  #ref     ACTG  ACTG   ACTG
  m_alt = 'ACTGATAC' + 'ACTG'
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 1)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons.tolist(), [(3, 3, 2), (5, 7, -2)], beacons)

  chrom = [(0, 0), (1, 1)]
  #        12345678   9012
  #ref     ACTGACTG   ACTG
  m_alt = 'ACTGAC' + 'ACTG'
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 1)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons.tolist(), [(5, 5, -2)], beacons)

  pos = [2, 6, 8]
  stop = [5, 7, 9]
  ref = ['CTG', 'C', 'G']
  alt = ['C', 'T', 'GAT']
  p = [0.1, 0.1, 0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2), (1, 0), (2, 2)]
  #        12345678  9012
  #ref     ACTGACTG  ACTG
  m_alt = 'AC''ATTGATACTG'
  #        12  3456789012
  #         ^   ^ ^
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons.tolist(), [(1, 1, -2), (5, 3, 0), (7, 5, 2)])

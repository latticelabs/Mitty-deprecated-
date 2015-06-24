import numpy as np

import mitty.lib.variants as vr
import mitty.lib.reads as reads
from nose.tools import assert_sequence_equal


def expand_seq_test1():
  """Expand sequence: no variants"""
  ref_seq = 'ACTGACTGACTGACT'
  ml = vr.VariantList()
  chrom = []
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)

  assert ref_seq == alt_seq
  assert_sequence_equal(beacons[1:-1], [])


def expand_seq_test2():
  """Expand sequence: one variant of each type"""
  #          012345678901234
  ref_seq = 'ACTGACTGACTGACT'

  pos = [3]
  stop = [4]
  ref = ['G']
  alt = ['T']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 1)]
  #        012345678901234
  #ref     ACTGACTGACTGACT
  m_alt = 'ACTTACTGACTGACT'
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert ref_seq == alt_seq, alt_seq
  assert_sequence_equal(beacons[1:-1], [], str(beacons))

  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 1)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(3, 3, 0)])

  pos = [3]
  stop = [4]
  ref = ['G']
  alt = ['GAAA']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 0)]

  #        0123   45678901234
  #ref     ACTG   ACTGACTGACT
  m_alt = 'ACTGAAAACTGACTGACT'
  #        012345678901234567
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(4, 4, 3)])

  pos = [3]
  stop = [7]
  ref = ['GACT']
  alt = ['G']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2)]

  #        012345678901234
  #            xxx
  #ref     ACTGACTGACTGACT
  m_alt = 'ACTG' 'GACTGACT'
  #        0123   45678901
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(7, 4, -3)])


def expand_seq_test3():
  """Expand sequence: multiple variants, different combinations"""
  #          012345678901
  ref_seq = 'ACTGACTGACTG'

  pos = [3, 5]
  stop = [4, 8]
  ref = ['G', 'CTG']
  alt = ['GAT', 'C']
  p = [0.1, 0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 1), (1, 1)]
  #        0123  4567   8901
  #ref     ACTG  ACTG   ACTG
  m_alt = 'ACTGATAC' + 'ACTG'
  #        01234567     8901
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 1)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(4, 4, 2), (8, 8, -2)])

  chrom = [(0, 0), (1, 1)]
  #        012345678901
  #ref     ACTGACTGACTG
  m_alt = 'ACTGAC''ACTG'
  #        012345  6789
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 1)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(8, 6, -2)])

  pos = [1, 5, 7]
  stop = [4, 6, 8]
  ref = ['CTG', 'C', 'G']
  alt = ['C', 'T', 'GAT']
  p = [0.1, 0.1, 0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2), (1, 0), (2, 2)]
  #        01234567  8901
  #ref     ACTGACTG  ACTG
  m_alt = 'AC''ATTGATACTG'
  #        01  2345678901
  #         ^   ^ ^
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(4, 2, -2), (5, 3, 0), (8, 6, 2)])


def expand_seq_test4():
  """Expand sequence: different variants at same POS"""
  #          012345
  ref_seq = 'AAAGGG'

  pos = [2, 2]
  stop = [3, 3]
  ref = ['A', 'A']
  alt = ['ATT', 'C']
  p = [0.1, 0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 0), (1, 1)]
  #        012  345
  #ref     AAA  GGG
  m_alt = 'AAATTGGG'  # For copy 0
  #        01234567
  alt_seq, beacons = reads.expand_sequence(ref_seq, ml, chrom, 0)
  assert alt_seq == m_alt, alt_seq
  assert_sequence_equal(beacons[1:-1].tolist(), [(3, 3, 2)])


def cigar_test1():
  """Rolling cigars: No variants"""
  #          012345678901234
  ref_seq = 'ACTGACTGACTGACT'
  ml = vr.VariantList()
  chrom = []
  alt_seq, variant_waypoints = reads.expand_sequence(ref_seq, ml, chrom, 0)
  read_list = np.rec.fromrecords([(n, 3) for n in range(0, 8)], names=['start_a', 'read_len'])
  pos, cigars = reads.roll_cigars(variant_waypoints, read_list)
  for cigar in cigars:
    assert cigar == '3=', cigar
  for n, p in enumerate(pos):
    assert p == n


def cigar_test2():
  """Rolling cigars: SNP"""
  #          012345678901234
  ref_seq = 'ACTGACTGACTGACT'
  #             ^
  pos = [3]
  stop = [4]
  ref = ['G']
  alt = ['T']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2)]

  #      012345678901234
  #ref   ACTGACTGACTGACT
  #alt   ACTTACTGACTGACT
  #         ^
  alt_seq, variant_waypoints = reads.expand_sequence(ref_seq, ml, chrom, 0)
  read_list = np.rec.fromrecords([(n, 3) for n in range(6)], names=['start_a', 'read_len'])
  pos, cigars = reads.roll_cigars(variant_waypoints, read_list)
  assert cigars[0] == '3=', cigars[0]
  assert cigars[1] == '2=1X', cigars[1]
  assert cigars[2] == '1=1X1=', cigars[2]
  assert cigars[3] == '1X2=', cigars[3]
  assert cigars[4] == '3=', cigars[4]
  assert cigars[5] == '3=', cigars[5]

  for n, p in enumerate(pos):
    assert p == n


def cigar_test3():
  """Rolling cigars: INS (incl. soft-clip)"""
  ref_seq = 'ACTGACTGACTGACT'
  pos = [3]
  stop = [4]
  ref = ['G']
  alt = ['GA']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2)]

  #        0123 45678901234
  #ref     ACTG ACTGACTGACT
  #alt     ACTGAACTGACTGACT
  #        0123456789012345
  alt_seq, variant_waypoints = reads.expand_sequence(ref_seq, ml, chrom, 0)
  read_list = np.rec.fromrecords([(n, 3) for n in range(6)], names=['start_a', 'read_len'])
  pos, cigars = reads.roll_cigars(variant_waypoints, read_list)
  assert cigars[0] == '3=', cigars[0]
  assert cigars[1] == '3=', cigars[1]
  assert cigars[2] == '2=1S', cigars[2]
  assert cigars[3] == '1=1I1=', cigars[3]
  assert cigars[4] == '1S2=', cigars[4]
  assert cigars[5] == '3=', cigars[5]

  assert pos[0] == 0
  assert pos[1] == 1
  assert pos[2] == 2
  assert pos[3] == 3
  assert pos[4] == 4
  assert pos[5] == 4


def cigar_test4():
  """Rolling cigars: read completely inside INS"""
  ref_seq = 'ACTGACTGACTGACT'
  pos = [3]
  stop = [4]
  ref = ['G']
  alt = ['GAAAA']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2)]

  #        0123    45678901234
  #ref     ACTG    ACTGACTGACT
  #alt     ACTGAAAAACTGACTGACT
  #        0123456789012345678
  alt_seq, variant_waypoints = reads.expand_sequence(ref_seq, ml, chrom, 0)
  read_list = np.rec.fromrecords([(n, 3) for n in range(7)], names=['start_a', 'read_len'])
  pos, cigars = reads.roll_cigars(variant_waypoints, read_list)
  assert cigars[0] == '3=', cigars[0]
  assert cigars[1] == '3=', cigars[1]
  assert cigars[2] == '2=1S', cigars[2]
  assert cigars[3] == '1=2S', cigars[3]
  assert cigars[4] == '3S', cigars[4]
  assert cigars[5] == '3S', cigars[5]
  assert cigars[6] == '2S1=', cigars[6]

  assert pos[0] == 0
  assert pos[3] == 3
  assert pos[4] == 4
  assert pos[6] == 4


def cigar_test5():
  """Rolling cigars: DEL (incl. start at deletion)"""
  ref_seq = 'ACTGACTGACTGACT'
  pos = [3]
  stop = [6]
  ref = ['GAC']
  alt = ['G']
  p = [0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2)]

  #        012345678901234
  #ref     ACTGACTGACTGACT
  #alt     ACTG  TGACTGACT
  #        0123  456789012
  alt_seq, variant_waypoints = reads.expand_sequence(ref_seq, ml, chrom, 0)
  read_list = np.rec.fromrecords([(n, 3) for n in range(6)], names=['start_a', 'read_len'])
  pos, cigars = reads.roll_cigars(variant_waypoints, read_list)
  assert cigars[0] == '3=', cigars[0]
  assert cigars[1] == '3=', cigars[1]
  assert cigars[2] == '2=2D1=', cigars[2]
  assert cigars[3] == '1=2D2=', cigars[3]
  assert cigars[4] == '3=', cigars[4]
  assert cigars[5] == '3=', cigars[5]

  assert pos[0] == 0
  assert pos[3] == 3
  assert pos[4] == 6, pos
  assert pos[5] == 7


def cigar_test6():
  """Rolling cigars: SNP, INS, DEL in different combinations"""
  #          012345678901234
  ref_seq = 'ACTGACTGACTGACT'

  #       0  1234567890  1234
  #ref    A  CTGACTGACT  GACT
  #alt    AAACGGA  GTCTTTGACT
  #       0123456  7890123456
  #           ^     ^
  pos = [0, 2, 4, 8, 10]
  stop = [1, 3, 7, 9, 11]
  ref = ['A', 'T', 'ACT', 'A', 'T']
  alt = ['AAA', 'G', 'A', 'T', 'TTT']
  p = [0.1, .1, .1, .1, .1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 2), (1, 2), (2, 2), (3, 2), (4, 2)]

  alt_seq, variant_waypoints = reads.expand_sequence(ref_seq, ml, chrom, 0)
  read_list = np.rec.fromrecords([(0, 3), (0, 4), (0, 5), (0, 7), (1, 7), (4, 8)],
                                 names=['start_a', 'read_len'])
  pos, cigars = reads.roll_cigars(variant_waypoints, read_list)
  assert cigars[0] == '1=2S', cigars[0]
  assert cigars[1] == '1=2I1=', cigars[1]
  assert cigars[2] == '1=2I1=1X', cigars[2]
  assert cigars[3] == '1=2I1=1X2=', cigars[3]
  assert cigars[4] == '2S1=1X2=2D1=', cigars[4]
  assert cigars[5] == '1X2=2D1=1X2=1S', cigars[5]


def cigar_test7():
  """Rolling cigars: different variants at same POS"""
  #          012345
  ref_seq = 'AAAGGG'

  pos = [2, 2]
  stop = [3, 3]
  ref = ['A', 'A']
  alt = ['ATT', 'C']
  p = [0.1, 0.1]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 0), (1, 1)]
  #        012  345
  #ref     AAA  GGG
  m_alt = 'AAATTGGG'  # For copy 0
  #        01234567
  alt_seq, variant_waypoints = reads.expand_sequence(ref_seq, ml, chrom, 0)
  read_list = np.rec.fromrecords([(0, 3), (1, 3), (2, 3), (3, 3), (4, 3)],
                                 names=['start_a', 'read_len'])
  pos, cigars = reads.roll_cigars(variant_waypoints, read_list)
  assert cigars[0] == '3=', cigars[0]
  assert cigars[1] == '2=1S', cigars[1]
  assert cigars[2] == '1=2S', cigars[2]
  assert cigars[3] == '2S1=', cigars[3]
  assert cigars[4] == '1S2=', cigars[4]

  #        012345
  #ref     AAAGGG
  m_alt = 'AAACGG'  # For copy 1
  #        012345

  alt_seq, variant_waypoints = reads.expand_sequence(ref_seq, ml, chrom, 1)
  read_list = np.rec.fromrecords([(0, 3), (1, 3), (2, 3), (3, 3)],
                                 names=['start_a', 'read_len'])
  pos, cigars = reads.roll_cigars(variant_waypoints, read_list)
  assert cigars[0] == '3=', cigars[0]
  assert cigars[1] == '3=', cigars[1]
  assert cigars[2] == '3=', cigars[2]
  assert cigars[3] == '3=', cigars[3]


def old_style_cigar_test():
  """Converting extended cigars to old style cigars"""
  assert reads.old_style_cigar('100=') == '100M'
  assert reads.old_style_cigar('20=1X40=') == '61M'
  assert reads.old_style_cigar('20S1X40=') == '20S41M'
  assert reads.old_style_cigar('20S1X30I40=') == '20S1M30I40M'
  assert reads.old_style_cigar('20S1X30I40=30D1X1X20M') == '20S1M30I40M30D22M'
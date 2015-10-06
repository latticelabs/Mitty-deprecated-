from nose.tools import assert_raises

import mitty.benchmarking.creed as creed


def categorized_reads_class_test():
  """benchmarking: CategorizedReads class"""
  assert_raises(RuntimeError, creed.CategorizedReads)  # Need seq len and seq id for new file
  seq_names = ['chr1', 'chr2']
  seq_lengths = [100, 200]
  cat_reads = creed.CategorizedReads(seq_names=seq_names, seq_lengths=seq_lengths, copies=2)
  assert cat_reads.get_sequence_info()[1][1] == seq_lengths[1]

  cat_reads.append(1, 0, 20, 100, 0b10010)
  cat_reads.append(1, 0, 10, 100, 0b10000)
  cat_reads.append(0, 1, 10, 100, 0b10000)
  cat_reads.finalize()

  dta = cat_reads.get_data(1, 0)
  assert dta.shape[0] == 2
  assert dta['pos'][0] < dta['pos'][1]

  assert_raises(RuntimeError, cat_reads.append, 1, 0, 10, 100, 0b10000)  # Can't append after finalizing


class MyRead:
  def __init__(self, qname, secondary, paired, read1, unmapped, reference_id, pos, cigarstring):
    self.qname, self.is_secondary, self.is_paired, self.is_read1 = qname, secondary, paired, read1
    self.is_unmapped, self.reference_id, self.pos, self.cigarstring = unmapped, reference_id, pos, cigarstring


def analyze_read_test():
  """Read analysis"""
  qname = '3|15|0|1|898|100|100=|0|744|100|24=2I74='
  read = MyRead(qname=qname, secondary=False, paired=True, read1=True, unmapped=False, reference_id=14, pos=898, cigarstring='100M')
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped = creed.analyze_read(read, window=0, extended=False)
  assert read_serial == 30, read_serial
  assert chrom == 15, chrom
  assert cpy == 0, cpy
  assert ro == 1, ro
  assert pos == 898, pos
  assert cigar == '100M', cigar
  assert chrom_c and pos_c and cigar_c == 1

  read = MyRead(qname=qname, secondary=False, paired=True, read1=False, unmapped=False, reference_id=14, pos=744, cigarstring='24M2I74M')
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped = creed.analyze_read(read, window=0, extended=False)
  assert read_serial == 31, read_serial
  assert chrom == 15, chrom
  assert cpy == 0, cpy
  assert ro == 0, ro
  assert pos == 744, pos
  assert cigar == '24M2I74M', cigar
  assert chrom_c and pos_c and cigar_c == 1

  read = MyRead(qname=qname, secondary=False, paired=True, read1=True, unmapped=False, reference_id=14, pos=898, cigarstring='24M2I74M')
  read_serial, chrom, cpy, ro, pos, rl, cigar, ro_m, pos_m, rl_m, cigar_m, chrom_c, pos_c, cigar_c, unmapped = creed.analyze_read(read, window=0, extended=False)
  assert read_serial == 30, read_serial
  assert chrom == 15, chrom
  assert cpy == 0, cpy
  assert ro == 1, ro
  assert pos == 898, pos
  assert cigar == '100M', cigar
  assert chrom_c == 1
  assert pos_c == 1
  assert cigar_c == 0

  # # Test if we can do fuzzy matching if we soft clips and so on
  # qname = '3|15|0|0|100|1M1000D99M|1|200|100='
  # read = MyRead(qname=qname, secondary=False, paired=True, read1=True, unmapped=False, reference_id=14, pos=180, cigarstring='1S99M')
  # read_serial, chrom, cpy, ro, pos, cigar, read_category = creed.analyze_read(read, window=0, extended=False)
  # assert read_category == 0b100100, bin(read_category)


# def create_sample_misaligned_bam(per_bam_fname):
#   """Utility function to create a sample BAM with specific read errors. Return us the feature positions and the
#   correct answers for the three types of read categories as creed.count_reads_under_features should return from
#   this file and the feature positions. See sketch in test folder"""
#
#   feature_positions = [  # 0,1 = het 2 = hom
#     (10, 20, 0), (70, 80, 2), (90, 91, 2), (110, 120, 1), (140, 141, 2), (150, 160, 1), (170, 180, 2), (210, 240, 0)]
#
#   # ((correct start, correct stop, aligned start, aligned stop),
#   #  ( ... mate ...)), ...
#
#   read_positions = [
#     [  # For copy 0
#       ((15, 30, 15, 30),
#        (40, 50, 45, 55)),
#       ((55, 65, 55, 65),
#        (75, 85, 75, 85)),
#       ((75, 95, 75, 95),
#        (100, 110, 200, 220)),
#       ((100, 110, 150, 160),
#        (120, 130, 120, 130)),
#       ((135, 150, 135, 150)
#        (160, 175, 160, 175)),
#       ((145, 155, 145, 155),
#        (165, 175, 165, 175)),
#       ((195, 200, 195, 200),
#        (215, 225, 230, 240),
#        (230, 235, 230, 235))],  # Copy 0  ...
#     [  # Copy 1
#       ((15, 30, 15, 30),
#        (40, 50, 40, 50)),
#       ((55, 65, 55, 65),
#        (75, 80, 75, 80)),
#       ((85, 95, 85, 95),
#        (105, 115, 105, 115)),
#       ((152, 158, 152, 158),
#        (172, 178, 172, 178))]  # Copy 1
#   ]
#
#   import pysam
#   bam_fp = pysam.AlignmentFile(per_bam_fname, 'wb')
#   for r_pos in read_positions:
#     r = pysam.AlignedSegment()



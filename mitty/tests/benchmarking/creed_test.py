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
  qname = '3|15|0|1|898|100=|0|744|24=2I74='
  read = MyRead(qname=qname, secondary=False, paired=True, read1=True, unmapped=False, reference_id=14, pos=898, cigarstring='100M')
  read_serial, chrom, cpy, ro, pos, cigar, read_category = creed.analyze_read(read, window=0, extended=False)
  assert read_serial == 30, read_serial
  assert chrom == 15, chrom
  assert cpy == 0, cpy
  assert ro == 1, ro
  assert pos == 898, pos
  assert cigar == '100M', cigar
  assert read_category == 0b010000, bin(read_category)

  read = MyRead(qname=qname, secondary=False, paired=True, read1=False, unmapped=False, reference_id=14, pos=744, cigarstring='24M2I74M')
  read_serial, chrom, cpy, ro, pos, cigar, read_category = creed.analyze_read(read, window=0, extended=False)
  assert read_serial == 31, read_serial
  assert chrom == 15, chrom
  assert cpy == 0, cpy
  assert ro == 0, ro
  assert pos == 744, pos
  assert cigar == '24M2I74M', cigar
  assert read_category == 0b100000, bin(read_category)

  read = MyRead(qname=qname, secondary=False, paired=True, read1=True, unmapped=False, reference_id=14, pos=898, cigarstring='24M2I74M')
  read_serial, chrom, cpy, ro, pos, cigar, read_category = creed.analyze_read(read, window=0, extended=False)
  assert read_serial == 30, read_serial
  assert chrom == 15, chrom
  assert cpy == 0, cpy
  assert ro == 1, ro
  assert pos == 898, pos
  assert cigar == '100M', cigar
  assert read_category == 0b010100, bin(read_category)

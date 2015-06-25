import pysam
from nose.tools import assert_raises

import mitty.benchmarking.perfectbam as pbam


def categorized_reads_class_test():
  """benchmarking: CategorizedReads class"""
  assert_raises(RuntimeError, pbam.CategorizedReads)  # Need seq len and seq id for new file
  seq_names = ['chr1', 'chr2']
  seq_lengths = [100, 200]
  cat_reads = pbam.CategorizedReads(seq_names=seq_names, seq_lengths=seq_lengths, copies=2)
  assert cat_reads.get_sequence_info()[1][1] == seq_lengths[1]

  cat_reads.append(1, 0, 10, 100, 0b10000)
  cat_reads.append(1, 0, 20, 100, 0b10010)
  cat_reads.append(0, 1, 10, 100, 0b10000)
  cat_reads.finalize()

  dta = cat_reads.get_data(1, 0)
  assert dta.shape[0] == 2

  assert_raises(RuntimeError, cat_reads.append, 1, 0, 10, 100, 0b10000)  # Can't append after finalizing


def create_small_bam(bam_name):
  """Create a small bam file that has a correctly aligned read, a read with wrong chrom, a read with wrong pos and a
    read with wrong CIGAR. Make sure the analysis of this BAM is correct"""
  header = {
    'HD': {'SO': 'coordinate', 'VN': '1.3'},
    'PG': [{'CL': 'custom', 'ID': 'Mitty', 'PN': 'Mitty', 'VN': '1'}],
      'SQ': [
        {'LN': 10, 'SN': 'chr1'},
        {'LN': 20, 'SN': 'chr2'}
      ]
  }
  with pysam.AlignmentFile(bam_name, 'wb', header=header) as bam_out_fp:
    # Correctly aligned read
    r1 = pysam.AlignedSegment()
    r1.query_name = 'r0|1|0|0|1|4=|1|5|4='
    r1.query_sequence = 'ACTG'
    r1.flag = 0x1 | 0x2 | 0x20 | 0x40
    r1.reference_id = 0
    r1.pos = 0
    r1.cigar = '4M'

    # Wrong chrom, correct pos, correct CIGAR
    r2 = pysam.AlignedSegment()
    r2.query_name = 'r0|1|0|0|1|4=|1|5|4='
    r2.query_sequence = 'ACTG'
    r2.flag = 0x1 | 0x2 | 0x10 | 0x80
    r2.reference_id = 1
    r2.pos = 5
    r2.cigar = '4M'

    # Correct chrom, wrong pos
    r3 = pysam.AlignedSegment()
    r3.query_name = 'r0|1|0|0|1|4=|1|5|4='
    r3.query_sequence = 'ACTG'
    r3.flag = 0b000001100011
    r3.reference_id = 0
    r3.pos = 3
    r3.cigar = '4M'

    # Correct chrom, correct pos, wrong CIGAR
    r4 = pysam.AlignedSegment()
    r4.query_name = 'r0|1|0|0|1|4=|1|5|4='
    r4.query_sequence = 'ACTG'
    r4.flag = 0b000010010011
    r4.reference_id = 1
    r4.pos = 5
    r4.cigar = '4M'


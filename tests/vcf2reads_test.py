from tests import *
from mitty.vcf2reads import *
from numpy.testing import assert_equal

def generator_test1():
  """Sequence expand no variants."""
  seq = 'ACTGACTGACTGACT'
  vg = get_variant_sequence_generator(ref_chrom_seq=seq, c1=[], chrom_copy=0, block_len=4, over_lap_len=1)

  i, s, cs, a = vg.next()
  assert i == 0
  assert s == 'ACTG'
  assert cs == 'TGAC'

  i, s, cs, a = vg.next()
  assert i == 3
  assert s == 'GACTG'
  assert cs == 'CTGAC'

  i, s, cs, a = vg.next()
  assert i == 7
  assert s == 'GACTG'
  assert cs == 'CTGAC'

  i, s, cs, a = vg.next()
  assert i == 11
  assert s == 'GACT'
  assert cs == 'CTGA'


def generator_test2():
  """Sequence expand with one variant"""
  seq = 'ACTGACTGACTGACT'
  #     'ATGACTGACTGACT'
  c1 = [Variation(1, 3, 'AC', 'A', HOMOZYGOUS)]
  vg = get_variant_sequence_generator(ref_chrom_seq=seq, c1=c1, chrom_copy=0, block_len=4, over_lap_len=1)

  i, s, cs, a = vg.next()
  assert i == 0, i
  assert s == 'ATGA', s
  assert cs == 'TACT', cs
  assert_equal(a[0], [1, 3, 4, 5])
  assert_equal(a[1], [2, 1, 1, 1])

  i, s, cs, a = vg.next()
  assert i == 3, i
  assert s == 'ACTGA', s
  assert cs == 'TGACT', cs

  i, s, cs, a = vg.next()
  assert i == 7, i
  assert s == 'ACTGA', s
  assert cs == 'TGACT', cs

  i, s, cs, a = vg.next()
  assert i == 11, i
  assert s == 'ACT', s
  assert cs == 'TGA', cs
  assert_equal(a[0], [13, 14, 15])
  assert_equal(a[1], [1, 1, 1])


def generator_test3():
  """Sequence expand with two variants"""
  seq = 'ACTGACTGACTGACT'
  #     'ATGTTACTGACTGACT'
  c1 = [Variation(1, 3, 'AC', 'A', HET1),
        Variation(4, 5, 'G', 'GTT', HOMOZYGOUS)]
  vg = get_variant_sequence_generator(ref_chrom_seq=seq, c1=c1, chrom_copy=0, block_len=4, over_lap_len=1)

  i, s, cs, a = vg.next()
  assert i == 0, i
  assert s == 'ATGTT', s
  assert cs == 'TACAA', cs
  assert_equal(a[0], [1, 3, 4, 5, 5])
  assert_equal(a[1], [2, 1, 1, 0, 0])

  i, s, cs, a = vg.next()
  assert i == 4, i
  assert s == 'TACTG', s
  assert cs == 'ATGAC', cs
  assert_equal(a[0], [5, 5, 6, 7, 8])
  assert_equal(a[1], [0, 1, 1, 1, 1])

  i, s, cs, a = vg.next()
  assert i == 8, i
  assert s == 'GACTG', s
  assert cs == 'CTGAC', cs

  vg = get_variant_sequence_generator(ref_chrom_seq=seq, c1=c1, chrom_copy=1, block_len=4, over_lap_len=1)
  #     'ACTGTTACTGACTGACT'

  i, s, cs, a = vg.next()
  assert i == 0, i
  assert s == 'ACTGTT', s
  assert cs == 'TGACAA', cs


def reads_and_cigar_test():
  """Perfect read generator, no variants"""

  # ACTGGACTTGACCTGACTGAACTT
  # 012345678901234567890123
  #           1         2
  #     *          *
  #
  seq = 'ACTGGACTTGACCTGACTGAACTT'
  c_seq = seq.translate(DNA_complement)

  reads = perfect_reads_from_seq(ref=[seq, c_seq], ref_start=4, read_len=4, template_len=9, paired=True,
                                 strand_values=[[0, 0]])
  correct_reads = [
    [Read(4, '4M', 'GACT', '~~~~'), Read(9, '4M', 'GGTC', '~~~~')],
    [Read(4, '4M', 'GACT', '~~~~'), Read(9, '4M', 'GGTC', '~~~~')]
  ]
  assert reads == correct_reads, reads

  reads = perfect_reads_from_seq(ref=[seq, c_seq], ref_start=4, read_len=4, template_len=9, paired=True,
                                 strand_values=[[0, 0], [1], [], [0, 0, 0]])
  correct_reads = [
    [Read(4, '4M', 'GACT', '~~~~'), Read(9, '4M', 'GGTC', '~~~~')],  # x 2
    [Read(4, '4M', 'GACT', '~~~~'), Read(9, '4M', 'GGTC', '~~~~')],
    [Read(5, '4M', 'TGAA', '~~~~'), Read(10, '4M', 'TCCA', '~~~~')],  # x 1
    [Read(7, '4M', 'TTGA', '~~~~'), Read(12, '4M', 'TCAG', '~~~~')],  # x 3
    [Read(7, '4M', 'TTGA', '~~~~'), Read(12, '4M', 'TCAG', '~~~~')],
    [Read(7, '4M', 'TTGA', '~~~~'), Read(12, '4M', 'TCAG', '~~~~')]
  ]
  assert reads == correct_reads, reads


def perfect_reads_from_seq_test():
  """Perfect read generator, no variants"""

  # ACTGGACTTGACCTGACTGAACTT
  # 012345678901234567890123
  #           1         2
  #     *          *
  #
  seq = 'ACTGGACTTGACCTGACTGAACTT'
  c_seq = seq.translate(DNA_complement)

  reads = perfect_reads_from_seq(ref=[seq, c_seq], ref_start=4, read_len=4, template_len=9, paired=True,
                                 strand_values=[[0, 0]])
  correct_reads = [
    [Read(4, '4M', 'GACT', '~~~~'), Read(9, '4M', 'GGTC', '~~~~')],
    [Read(4, '4M', 'GACT', '~~~~'), Read(9, '4M', 'GGTC', '~~~~')]
  ]
  assert reads == correct_reads, reads

  reads = perfect_reads_from_seq(ref=[seq, c_seq], ref_start=4, read_len=4, template_len=9, paired=True,
                                 strand_values=[[0, 0], [1], [], [0, 0, 0]])
  correct_reads = [
    [Read(4, '4M', 'GACT', '~~~~'), Read(9, '4M', 'GGTC', '~~~~')],  # x 2
    [Read(4, '4M', 'GACT', '~~~~'), Read(9, '4M', 'GGTC', '~~~~')],
    [Read(5, '4M', 'TGAA', '~~~~'), Read(10, '4M', 'TCCA', '~~~~')],  # x 1
    [Read(7, '4M', 'TTGA', '~~~~'), Read(12, '4M', 'TCAG', '~~~~')],  # x 3
    [Read(7, '4M', 'TTGA', '~~~~'), Read(12, '4M', 'TCAG', '~~~~')],
    [Read(7, '4M', 'TTGA', '~~~~'), Read(12, '4M', 'TCAG', '~~~~')]
  ]
  assert reads == correct_reads, reads

generator_test2()
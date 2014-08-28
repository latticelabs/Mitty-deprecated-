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


def roll_cigar_test():
  """Rolling a lot of CIGARs"""
  # seq    ACTGA  CTG
  #r_seq   ACTGATTCTG
  #        1234566678
  pos_array = [[1, 2, 3, 4, 5, 6, 6, 6, 7, 8],
               [1, 1, 1, 1, 1, 0, 0, 1, 1, 1],
               [0, 0, 0, 0, 0, 1, 2, 0, 0, 0]]

  # Fully matching reads
  read_list = [[Read(_start_idx=0, _stop_idx=4)]]
  package_reads(read_list, pos_array)
  assert read_list[0][0].POS == 1
  assert read_list[0][0].CIGAR == '4M'

  read_list = [[Read(_start_idx=1, _stop_idx=5)]]
  package_reads(read_list, pos_array)
  assert read_list[0][0].POS == 2
  assert read_list[0][0].CIGAR == '4M'

  # Reads with an insertion
  read_list = [[Read(_start_idx=2, _stop_idx=6)]]
  package_reads(read_list, pos_array)
  assert read_list[0][0].POS == 3
  assert read_list[0][0].CIGAR == '3M1S'

  read_list = [[Read(_start_idx=3, _stop_idx=7)]]
  package_reads(read_list, pos_array)
  assert read_list[0][0].POS == 4
  assert read_list[0][0].CIGAR == '2M2S', read_list[0].CIGAR

  read_list = [[Read(_start_idx=4, _stop_idx=8)]]
  package_reads(read_list, pos_array)
  assert read_list[0][0].POS == 5
  assert read_list[0][0].CIGAR == '1M2I1M'

  # Reads from inside an insertion
  read_list = [[Read(_start_idx=5, _stop_idx=6)]]
  package_reads(read_list, pos_array)
  assert read_list[0][0].POS == 5
  assert read_list[0][0].CIGAR == '>1'

  # Deletion
  # seq    ACTG
  #r_seq    CTG
  #        1234
  pos_array = [[2, 3, 4],
               [1, 1, 1],
               [0, 0, 0]]

  # Read with deletion at very start
  read_list = [[Read(_start_idx=0, _stop_idx=3)]]
  package_reads(read_list, pos_array)
  assert read_list[0][0].POS == 2
  assert read_list[0][0].CIGAR == '3M'

  # seq    ACTG
  #r_seq   A TG
  #        1234
  pos_array = [[1, 3, 4],
               [2, 1, 1],
               [0, 0, 0]]

  # Read with deletion in middle
  read_list = [[Read(_start_idx=0, _stop_idx=3)]]
  package_reads(read_list, pos_array)
  assert read_list[0][0].POS == 1
  assert read_list[0][0].CIGAR == '1M1D2M', read_list[0][0].CIGAR

  # seq    ACTGA
  #r_seq   ACT A
  #        12345
  pos_array = [[1, 2, 3, 5],
               [1, 1, 2, 1],
               [0, 0, 0, 0]]

  # Read with deletion right at the end
  read_list = [[Read(_start_idx=0, _stop_idx=3)]]
  package_reads(read_list, pos_array)
  assert read_list[0][0].POS == 1
  assert read_list[0][0].CIGAR == '3M', read_list[0][0].CIGAR


# class Read_Module():
#   """A simple read model simulator for testing purposes."""
#   @staticmethod
#   def initialize(model_params, master_seed):
#     return {
#       'max_read_len': model_params['max_read_len'],
#       'read_advance': model_params['read_advance'],
#       'last_idx': 0
#     }
#
#   @staticmethod
#   def max_read_len(read_model_state):
#     return read_model_state['max_read_len']
#
#   @staticmethod
#   def overlap_len(read_model_state):
#     return read_model_state['max_read_len'] - 1
#
#
#   @staticmethod
#   def generate_reads(this_idx, this_seq_block, this_c_seq_block, this_arr, read_model_state):
#     template_list = []
#     lidx = read_model_state['last_idx']
#     ra = read_model_state['read_advance']
#     mrl = model_params['max_read_len']
#     for n in range(lidx + ra - this_idx, len(this_seq_block) - mrl, ra):
#       [Read(perfect_seq=this_seq_block[n:n+mrl], _start_idx=n, _stop_idx=n+mrl),
#        Read(perfect_seq=this_c_seq_block[n+mrl:n-1:-1], _start_idx=n, _stop_idx=n+mrl)]
#
#
#
#
#     return \
#       [
#       ]


class Read_Module():
  """A simple read model simulator for testing purposes."""
  @staticmethod
  def initialize(model_params, master_seed):
    return {
      'read_len': model_params['read_len'],
      'template_len': model_params['template_len'],
      'read_advance': model_params['read_advance'],
      'last_idx': -model_params['read_advance']
    }

  @staticmethod
  def max_read_len(read_model_state):
    return read_model_state['read_len']

  @staticmethod
  def overlap_len(read_model_state):
    return read_model_state['template_len'] - 1


  @staticmethod
  def generate_reads(this_idx, this_seq_block, this_c_seq_block, this_arr, read_model_state, master_seed):
    template_list = []
    l_idx = read_model_state['last_idx']
    ra = read_model_state['read_advance']
    rl = read_model_state['read_len']
    tl = read_model_state['template_len']
    for n in range(l_idx + ra - this_idx, len(this_seq_block) - tl, ra):
      template_list.append(
      [Read(perfect_seq=this_seq_block[n:n+rl], _start_idx=n, _stop_idx=n+rl),
       Read(perfect_seq=this_c_seq_block[n+tl-1:n+tl-rl-1:-1], _start_idx=n+tl-rl, _stop_idx=n+tl)]
      )
    read_model_state['last_idx'] = n
    return template_list, read_model_state


def reads_from_genome_test1():
  """Read generator basic test"""
  read_model = Read_Module()
  model_params = {'read_len': 2, 'template_len': 5, 'read_advance': 2}
  ref = {1: 'ACTGACTG'}
  chrom_list = [1, 2]  # The '2' should return us no reads, as genome has no chromosome 2
  # Generate variant sequence in one go
  rg = reads_from_genome(ref=ref, g1={}, chrom_list=chrom_list,
                         read_model=read_model, model_params=model_params,
                         block_len=10, master_seed=1)
  tl, chrom, cc = rg.next()
  assert tl[0][0].POS == 1, tl
  assert tl[0][1].perfect_seq == 'TC', tl

  assert tl[1][0].POS == 3, tl
  assert tl[1][1].perfect_seq == 'AG', tl

  # Sequence in blocks, results should be the same
  rg = reads_from_genome(ref=ref, g1={}, chrom_list=chrom_list,
                         read_model=read_model, model_params=model_params,
                         block_len=6, master_seed=1)
  tl, chrom, cc = rg.next()
  assert tl[0][0].POS == 1, tl
  assert tl[0][1].perfect_seq == 'TC', tl

  tl, chrom, cc = rg.next()
  assert tl[0][0].POS == 3, tl
  assert tl[0][1].perfect_seq == 'AG', tl


def write_test():
  import io
  fp, fp_c = io.BytesIO(), io.BytesIO()
  template_list = [  # Unpaired reads
    [Read(1, '2M', 'AC', 'AT', '~~')],
    [Read(10, '2M', 'AG', 'TG', '~!')]
  ]
  chrom = 1
  cc = 0
  serial_no = 1
  write_reads_to_file(fp, fp_c, template_list, chrom, cc, serial_no, write_corrupted=True)

  fp.seek(0)
  fastq = fp.read()
  correct_fastq = """@1:0|r1|1|2M
AC
+
~~
@1:0|r2|10|2M
AG
+
~~\n"""
  assert fastq == correct_fastq, fastq

  fp_c.seek(0)
  fastq = fp_c.read()
  correct_fastq = """@1:0|r1|1|2M
AT
+
~~
@1:0|r2|10|2M
TG
+
~!\n"""
  assert fastq == correct_fastq, fastq

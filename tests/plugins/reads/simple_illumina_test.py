from mitty.plugins.reads.simple_illumina_plugin import *


def test_initialize():
  """Model initialization"""
  model_params = {
    'paired': True,      # Are the reads paired or not.
    'read_len': 10,     # length of each read
    'template_len': 25, # length of template (only used for paired reads)
    'coverage': 1.0,     # How many x coverage do we want
    'max_p_error': 0.01, # Maximum error rate at tip of read
    'k': 0.3
  }
  model_state = initialize(model_params)
  assert 'coverage' in model_state
  assert model_state['error_profile'][9] == 0.01
  assert model_state['error_profile'][0] == 0.01 * 0.3 ** 9

  assert model_state['PHRED'][9] == chr(int(33 + -10*numpy.log10(0.01)))


def test_extract_reads_se():
  """Perfect read generator se"""
  seq =   'ACTGACTGACTGACTG'
  seq_c = 'TGACTGACTGACTGAC'
  t_start = numpy.array([0, 2, 4, 6])
  strand =  numpy.array([0, 1, 0, 1])
  rl = 5
  tl = 11
  template_list = extract_reads([seq, seq_c], t_start, strand, tl, rl, paired=False)

  assert len(template_list) == 4
  assert len(template_list[0]) == 1
  assert template_list[0][0].perfect_seq == 'ACTGA'
  assert template_list[1][0].perfect_seq == 'AGTCA'


def test_extract_reads_paired():
  """Perfect read generator paired"""
  seq =   'ACTGACTGACTGACTG'
  seq_c = 'TGACTGACTGACTGAC'
  t_start = numpy.array([0, 2, 4, 6])
  strand =  numpy.array([0, 1, 0, 1])
  rl = 5
  tl = 11
  template_list = extract_reads([seq, seq_c], t_start, strand, tl, rl, paired=True)

  assert len(template_list) == 4
  assert len(template_list[0]) == 2
  assert template_list[0][0].perfect_seq == 'ACTGA'
  assert template_list[0][1].perfect_seq == 'AGTCA'


def test_read_corruption_se():
  """Read corruption SE"""
  seq =   'ACTGACTGACTGACTG'
  seq_c = 'TGACTGACTGACTGAC'
  t_start = numpy.array([0, 2, 4, 6])
  strand =  numpy.array([0, 1, 0, 1])
  rl = 5
  tl = 11
  idx = (numpy.array([0, 1, 1]), numpy.array([0, 0, 3]))
  corrupt_bases = 'GTA'
  template_list = extract_reads([seq, seq_c], t_start, strand, tl, rl, paired=False)
  fill_out_corrupt_bases(template_list, corrupt_bases, idx, '~~~~~')

  assert template_list[0][0].perfect_seq == 'ACTGA'
  assert template_list[0][0].corrupt_seq == 'GCTGA'

  assert template_list[1][0].perfect_seq == 'AGTCA'
  assert template_list[1][0].corrupt_seq == 'TGTAA'  # Note error locations


def test_read_corruption_paired():
  """Read corruption paired"""
  seq =   'ACTGACTGACTGACTG'
  seq_c = 'TGACTGACTGACTGAC'
  t_start = numpy.array([0, 2, 4, 6])
  strand =  numpy.array([0, 1, 0, 1])
  rl = 5
  tl = 11
  idx = (numpy.array([0, 1, 1]), numpy.array([0, 0, 3]))
  corrupt_bases = 'GTA'
  template_list = extract_reads([seq, seq_c], t_start, strand, tl, rl, paired=True)
  fill_out_corrupt_bases(template_list, corrupt_bases, idx)

  assert template_list[0][0].perfect_seq == 'ACTGA'
  assert template_list[0][0].corrupt_seq == 'GCTGA'

  assert template_list[0][1].perfect_seq == 'AGTCA'
  assert template_list[0][1].corrupt_seq == 'TGTAA'  # Note error locations

"""This is the stock read plugin that approximates illumina reads

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
__example_param_text = """
{
  'read_len': 100,     #length of each read
  'template_len': 250, #length of template
  'coverage': 1.0,     #How many x coverage do we want
  'max_p_error': 0.01,  #Maximum error rate at tip of read
  'k': 0.3
}
"""

_description = """
This read generator generates Illumina like reads with a exponential error profile
Example parameter set:
""" + __example_param_text

_example_params = eval(__example_param_text)

from mitty.lib.read import Read, direction
import numpy
import logging
logger = logging.getLogger(__name__)


def initialize(model_params):
  """The only interesting this we do here is setup the error profile."""
  state = dict(model_params)
  state['error_profile'] = [state['max_p_error'] * state['k'] ** n for n in range(state['read_len'])][::-1]
  state['PHRED'] = ''.join([chr(int(33 + max(0, min(-10*numpy.log10(p), 93)))) for p in state['error_profile']])
  return state


def max_read_len(read_model_state):
  return read_model_state['read_len']


def overlap_len(read_model_state):
  return read_model_state['template_len'] - 1


def init_rngs(master_seed):
  read_loc_rng_seed, read_order_rng_seed, error_loc_rng_seed, base_choice_rng_seed = \
    numpy.random.RandomState(seed=master_seed).randint(100000000, size=4)
  logger.debug('Used master seed to generate seeds {:d}, {:d}, {:d}, {:d}'.
               format(read_loc_rng_seed, read_order_rng_seed, error_loc_rng_seed, base_choice_rng_seed))
  return {
    'loc_rng': numpy.random.RandomState(seed=read_loc_rng_seed).randint,
    'order_rng': numpy.random.RandomState(seed=read_order_rng_seed).randint,
    'error_loc_rng': numpy.random.RandomState(seed=error_loc_rng_seed).rand,
    'base_choice_rng': numpy.random.RandomState(base_choice_rng_seed).choice
  }


def generate_reads(this_idx, this_seq_block, this_c_seq_block, this_arr, read_model_state, master_seed):
  rngs = init_rngs(master_seed)
  seqs = [this_seq_block, this_c_seq_block]  # Just package it nicely
  template_count = 0.5 * read_model_state['coverage'] * len(this_seq_block) / (read_model_state['read_len'] * 2.0)
  rl = read_model_state['read_len']
  tl = read_model_state['template_len']
  t_start = rngs['loc_rng'](low=0, high=len(this_seq_block) - tl, size=template_count)
  read_order = rngs['order_rng'](2, size=template_count)

  template_list = extract_reads(seqs, t_start, read_order, tl, rl)

  if read_model_state['generate_corrupted_reads']:  # Corrupt the bases and fill out the corrupted_seq field
    erp = read_model_state['error_profile']
    rl = read_model_state['read_len']
    tc = len(template_list)
    idx = numpy.where(rngs['error_loc_rng'](tc * len(template_list[0]), rl) < erp)
    corrupt_bases = rngs['base_choice_rng'](['A','C','G','T'], size=idx[0].size, replace=True, p=[.3, .2, .2, .3]).tostring()

    fill_out_corrupt_bases(template_list, corrupt_bases, idx, read_model_state['PHRED'])

  return template_list


def extract_reads(seqs, t_start, read_order, tl, rl):
  """Refactored out random variables to make testing easier."""
  template_list = []
  for n in range(t_start.size):
    r1 = Read(perfect_seq=seqs[0][t_start[n]:t_start[n]+rl], _start_idx=t_start[n], _stop_idx=t_start[n]+rl, direction='>')
    r2 = Read(perfect_seq=seqs[1][t_start[n]+tl-rl:t_start[n]+tl][::-1], _start_idx=t_start[n]+tl-rl, _stop_idx=t_start[n]+tl, direction='<')
    if read_order[n]:
      template_list += [[r2, r1]]
    else:
      template_list += [[r1, r2]]
  return template_list


def fill_out_corrupt_bases(template_list, corrupt_bases, idx, phred):
  """Refactored out random variables to make testing easier."""
  rl = len(template_list[0][0].perfect_seq)
  r_phred = phred[::-1]

  idx_ptr = 0
  read_ctr = 0
  for template in template_list:
    for read in template:
      read.PHRED = phred if read.direction == '>' else r_phred
      if (idx_ptr == idx[0].size or  # We are out of read corruptions
          idx[0][idx_ptr] < read_ctr):  # or this read is not corrupted
        read.corrupt_seq = read.perfect_seq
      else:
        corrupt_seq = bytearray(read.perfect_seq)
        while idx[0][idx_ptr] == read_ctr:
          corrupt_seq[idx[1][idx_ptr]] = corrupt_bases[idx_ptr]  # Sequence ends always correspond to inner
          idx_ptr += 1
          if idx_ptr == idx[0].size:
            break
        read.corrupt_seq = corrupt_seq.__str__()
      read_ctr += 1


if __name__ == "__main__":
  print _description
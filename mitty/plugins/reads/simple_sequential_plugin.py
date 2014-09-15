"""Very simple deterministic read generator with no read corruption.

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
from mitty.lib.read import Read, direction

__example_param_text = """
{
  'paired': True,      #Are the reads paired or not.
  'read_len': 100,     #length of each read
  'template_len': 250, #length of template (only used for paired reads)
  'read_advance': 20   #how much to advance along the reference after generating a read. Determines "coverage"
}
"""

_description = """
This is a simple, deterministic read generator with no read corruption.
Example parameter set:
""" + __example_param_text

_example_params = eval(__example_param_text)


def initialize(model_params):
  return {
    'paired': model_params['paired'],
    'read_len': model_params['read_len'],
    'template_len': model_params['template_len'],
    'read_advance': model_params['read_advance']
  }


def max_read_len(read_model_state):
  return read_model_state['read_len']


def overlap_len(read_model_state):
  if read_model_state['paired']:
    return read_model_state['template_len'] - 1
  else:
    return read_model_state['read_len']


def generate_reads(this_idx, this_seq_block, this_c_seq_block, this_arr, read_model_data, master_seed):
  seq = [this_seq_block, this_c_seq_block]
  paired = read_model_data.get('paired', False)
  ra = read_model_data.get('read_advance', 10)
  rl = read_model_data.get('read_len', 100)
  tl = read_model_data.get('template_len', 250) if paired else rl
  template_list = []
  strand = 0
  for n in range(0, len(this_seq_block) - tl, ra):
    # n is in local coordinates - references to this_seq_block, this_arr etc.
    r1 = Read(perfect_seq=seq[1][n:n+rl][::-1] if strand else seq[0][n:n+rl],
              _start_idx=n, _stop_idx=n+rl, direction=direction[strand])
    if paired:
      r2 = Read(perfect_seq=seq[0][n+tl-rl:n+tl] if strand else seq[1][n+tl-rl:n+tl][::-1],
                _start_idx=n+tl-rl, _stop_idx=n+tl, direction=direction[1 - strand])
      template_list += [[r1, r2]]
    else:
      template_list += [[r1]]
    strand = 1 - strand
  return template_list

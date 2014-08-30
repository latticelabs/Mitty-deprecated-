from mitty.vcf2reads import Read


def initialize(model_params, master_seed):
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


def generate_reads(this_idx, this_seq_block, this_c_seq_block, this_arr, read_model_state, master_seed):
  paired = read_model_state['paired']
  ra = read_model_state['read_advance']
  rl = read_model_state['read_len']
  tl = read_model_state['template_len'] if paired else rl
  template_list = []
  for n in range(0, len(this_seq_block) - tl, ra):
    # n is in local coordinates - references to this_seq_block, this_arr etc.
    this_template = [Read(perfect_seq=this_seq_block[n:n+rl], _start_idx=n, _stop_idx=n+rl)]
    if paired:
     this_template += [Read(perfect_seq=this_c_seq_block[n+tl-1:n+tl-rl-1:-1], _start_idx=n+tl-rl, _stop_idx=n+tl)]
    template_list.append(this_template)
  return template_list, read_model_state

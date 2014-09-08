"""This is the stock read plugin that approximates illumina reads

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
__example_param_text = """
{
  'paired': True,      #Are the reads paired or not.
  'read_len': 100,     #length of each read
  'template_len': 250, #length of template (only used for paired reads)
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
  if read_model_state['paired']:
    return read_model_state['template_len'] - 1
  else:
    return read_model_state['read_len']


def init_rngs(master_seed):
  read_loc_rng_seed, read_strand_rng_seed, error_loc_rng_seed, base_choice_rng_seed = \
    numpy.random.RandomState(seed=master_seed).randint(100000000, size=4)
  logger.debug('Used master seed to generate seeds {:d}, {:d}, {:d}, {:d}'.
               format(read_loc_rng_seed, read_strand_rng_seed, error_loc_rng_seed, base_choice_rng_seed))
  return {
    'loc_rng': numpy.random.RandomState(seed=read_loc_rng_seed).randint,
    'strand_rng': numpy.random.RandomState(seed=read_strand_rng_seed).randint,
    'error_loc_rng': numpy.random.RandomState(seed=error_loc_rng_seed).rand,
    'base_choice_rng': numpy.random.RandomState(base_choice_rng_seed).choice
  }


def generate_reads(this_idx, this_seq_block, this_c_seq_block, this_arr, read_model_state, master_seed):
  rngs = init_rngs(master_seed)
  seqs = [this_seq_block, this_c_seq_block]  # Just package it nicely
  paired = read_model_state['paired']
  template_count = read_model_state['coverage'] * len(this_seq_block) / (read_model_state['read_len'] * 2.0 if paired else 1.0)
  rl = read_model_state['read_len']
  tl = read_model_state['template_len'] if paired else rl
  t_start = rngs['loc_rng'](low=0, high=len(this_seq_block) - tl, size=template_count)
  strand = rngs['strand_rng'](2, size=template_count)

  template_list = extract_reads(seqs, t_start, strand, tl, rl, paired)

  if read_model_state['generate_corrupted_reads']:  # Corrupt the bases and fill out the corrupted_seq field
    erp = read_model_state['error_profile']
    rl = read_model_state['read_len']
    tc = len(template_list)
    idx = numpy.where(rngs['error_loc_rng'](tc * len(template_list[0]), rl) < erp)
    corrupt_bases = rngs['base_choice_rng'](['A','C','G','T'], size=idx[0].size, replace=True, p=[.3, .2, .2, .3]).tostring()

    fill_out_corrupt_bases(template_list, corrupt_bases, idx, read_model_state['PHRED'])

  return template_list


def extract_reads(seqs, t_start, strand, tl, rl, paired):
  """Refactored out random variables to make testing easier."""
  template_list = []
  for n in range(t_start.size):
    r1 = Read(perfect_seq=seqs[1][t_start[n]:t_start[n]+rl][::-1] if strand[n] else seqs[0][t_start[n]:t_start[n]+rl],
              _start_idx=t_start[n], _stop_idx=t_start[n]+rl, direction=direction[strand[n]])
    if paired:
      r2 = Read(perfect_seq=seqs[0][t_start[n]+tl-rl:t_start[n]+tl] if strand[n] else seqs[1][t_start[n]+tl-rl:t_start[n]+tl][::-1],
                _start_idx=t_start[n]+tl-rl, _stop_idx=t_start[n]+tl, direction=direction[strand[n]])
      template_list += [[r1, r2]]
    else:
      template_list += [[r1]]
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
          if read.direction == '>':
            corrupt_seq[idx[1][idx_ptr]] = corrupt_bases[idx_ptr] # Sequence ends correspond to inner
          else:
            corrupt_seq[rl - 1 - idx[1][idx_ptr]] = corrupt_bases[idx_ptr]  # Sequence starts correspond to inner
          idx_ptr += 1
          if idx_ptr == idx[0].size:
            break
        read.corrupt_seq = corrupt_seq.__str__()
      read_ctr += 1

# def corrupt_reads(template_list, rngs, read_model_state):
#   erp = read_model_state['error_profile']
#   rl = read_model_state['read_len']
#   tc = len(template_list)
#   idx = numpy.where(rngs['error_loc_rng'](tc * len(template_list[0]), rl) > erp)
#   corrupt_bases = rngs['base_choice_rng'](['A','C','G','T'], size=idx[0].size, replace=True, p=[.3, .2, .2, .3])
#
#   idx_ptr = 0
#   read_ctr = 0
#   for template in template_list:
#     for read in template:
#       if (idx_ptr == idx[0].size or  # We are out of read corruptions
#           idx[0][idx_ptr] < read_ctr):  # or this read is not corrupted
#         read.corrupt_seq = read.perfect_seq
#       else:
#         corrupt_seq = bytearray(read.perfect_seq)
#         while idx[0][idx_ptr] == read_ctr:
#           corrupt_seq[idx[1][idx_ptr]] = corrupt_bases[idx_ptr]
#           idx_ptr += 1
#           if idx_ptr == idx[0].size:
#             break
#         read.corrupt_seq = corrupt_seq.__str__()
#       read_ctr += 1


# def average_read_len(read_len=None, **kwargs):
#   """Given the same parameters passed to generate_reads tell us what the average read len is going to be. reads.py
#   uses this in combination with coverage and seq_len to figure out how many reads we need."""
#   return read_len
#
#
# def max_read_len(read_len=None, **kwarg):
#   """Return maximum read length of reads."""
#   return read_len
#
#
# def read_generator(seq=None,
#                    read_start=0,
#                    read_stop=0,
#                    reads_per_call=1000,
#                    num_reads=10000,
#                    generate_corrupt_reads=False,
#                    paired=False,
#                    read_len=None,
#                    template_len=None,
#                    master_seed=None,
#                    read_loc_rng_seed=0,
#                    read_strand_rng_seed=1,
#                    error_rng_seed=1,
#                    base_chose_rng_seed=2,
#                    max_p_error=.8,
#                    k=.1,
#                    **kwargs):
#   """Given a sequence generate reads with the given characteristics. See note about RNG states. The RNGs are
#   initialized only on the first call to the generator. New calls to the generator will keep the same RNGs
#
#   Inputs:
#     seq              - pair of seq and complement_seq (stringlike) containing the DNA sequence to generate reads from
#     chrom_copy       - we use this to change the rng seeds to ensure we take different reads/different positions
#                        from each chromosome copy
#     read_start       - start generating reads from here (0.0, 1.0)
#     read_stop        - stop generating reads from here (0.0, 1.0)
#     reads_per_call   - each call to next will generate these many reads
#     num_reads        - how many reads do we want in total
#     paired           - paired reads or not
#     read_len         - Fixed read length (for this model)
#     template_len     - Fixed template length (for this model). Only needed if paired is True
#     read_loc_rng_seed- Seed for rng that drives the read location picker
#     max_p_error      - error probability for last base of read
#                        (0.0 is perfect reads, 1.0 -> every base is guaranteed to be wrong)
#     k                - exponential factor. 0 < k < 1 The closer this is to 0 the quicker the base error rate drops
#     error_loc_rng    - from generate_reads, contains the 2 RNGs we need
#     base_chose_rng   - we don't need to return them as they are passed by reference and the state is propagated
#                        back up to caller.
#     kwargs           - to swallow any other arguments
#
#   Outputs
#                                   _________ ( seq_str, quality_str, coordinate)
#     perfect_reads  -             /
#                  [
#                   [[ ... ], [ ... ]],
#                   [[ ... ], [ ... ]], -> inner list = 2 elements if paired reads, 1 otherwise
#                        .
#                        .
#                        .
#                  ] -> outer list = number of reads
#
#   Notes:
#   0. This yields a generator
#   1. Coordinate is 0-indexed
#   2. Quality: Sanger scale 33-126
#   3. The number of reads returned on each iteration can be less than reads_per_call because we toss out reads with Ns
#      in them
#   4. The only tricky thing is we store the rng state in a function attribute, so that with repeated calls to reads we
#      get fresh random numbers
#
#   Read corruption
#
#   1. We use a simple exponential model
#   2. We generate random numbers corresponding to base flips all at once to save time.
#
#   """
#   assert k <= 1.0
#   #'tis better to ask for forgiveness
#   try:
#     read_loc_rng_randint = read_generator.read_loc_rng_randint
#     read_strand_rng_randint = read_generator.read_strand_rng_randint
#     error_loc_rng_rand = read_generator.error_loc_rng_rand
#     base_chose_rng_choice = read_generator.base_chose_rng_choice
#     logger.debug('Using previously initialized RNGs')
#   except AttributeError:
#     # We need to initialize the rngs and a bunch of other stuff
#     if master_seed:
#       read_loc_rng_seed, read_strand_rng_seed, error_loc_rng_seed, base_choice_rng_seed = \
#         numpy.random.RandomState(seed=master_seed).randint(100000000, size=4)
#       logger.debug('Used master seed to generate seeds {:d}, {:d}, {:d}, {:d}'.
#                    format(read_loc_rng_seed, read_strand_rng_seed, error_loc_rng_seed, base_choice_rng_seed))
#
#     read_generator.read_loc_rng_randint = numpy.random.RandomState(seed=read_loc_rng_seed).randint
#     read_generator.read_strand_rng_randint = numpy.random.RandomState(seed=read_strand_rng_seed).randint
#     read_generator.error_loc_rng_rand = numpy.random.RandomState(seed=error_rng_seed).rand
#     read_generator.base_chose_rng_choice = numpy.random.RandomState(base_chose_rng_seed).choice
#
#     read_loc_rng_randint = read_generator.read_loc_rng_randint
#     read_strand_rng_randint = read_generator.read_strand_rng_randint
#     error_loc_rng_rand = read_generator.error_loc_rng_rand
#     base_chose_rng_choice = read_generator.base_chose_rng_choice
#
#   error_profile = [max_p_error * k ** n for n in range(read_len)][::-1]
#
#   rl = read_len
#   tl = template_len if paired else rl
#
#   if paired:
#     num_reads = max(1, num_reads / 2)
#     reads_per_call = max(1, reads_per_call / 2)
#
#   read_count = 0
#   if read_start + tl >= read_stop:  # We should raise StopIteration and quit
#     logger.error('Template len too large for given sequence.')
#     read_count = num_reads
#   while read_count < num_reads:
#     nominal_read_count = min(reads_per_call, num_reads - read_count)
#     rd_st = read_loc_rng_randint(low=read_start, high=read_stop - tl, size=nominal_read_count)
#     rd_strand = read_strand_rng_randint(2, size=nominal_read_count)
#     reads = []
#     if paired:
#       for this_rd_st, this_rd_stand in zip(rd_st, rd_strand):
#         # seq_1 = seq[this_rd_stand][this_rd_st:this_rd_st + rl]
#         # seq_2 = seq[1 - this_rd_stand][this_rd_st + tl - 1:this_rd_st + tl - 1 - rl:-1]
#         seq_1 = seq[0][this_rd_st:this_rd_st + rl]
#         seq_2 = seq[1][this_rd_st + tl - 1:this_rd_st + tl - 1 - rl:-1]
#         if 'N' in seq_1 or 'N' in seq_2:  # read taken from a masked/unknown region
#           continue
#         reads.append([[seq_1, '~' * rl, this_rd_st],
#                       [seq_2, '~' * rl, this_rd_st + tl - rl]])
#     else:
#       for this_rd_st in rd_st:
#         if 'N' in seq[this_rd_st:this_rd_st + tl]:  # read taken from a masked/unknown region
#           continue
#         reads.append([[seq[this_rd_st:this_rd_st + rl], '~' * rl, this_rd_st]])
#
#     if generate_corrupt_reads:
#       corrupted_reads = corrupt_reads(reads, error_profile, error_loc_rng_rand, base_chose_rng_choice)
#     else:
#       corrupted_reads = None
#
#     read_count += nominal_read_count * (2 if paired else 1) # len(reads) We don't use actual read count to avoid thrashing when we have tons of Ns in the sequence
#     logger.debug('Generated {:d} reads'.format(read_count))
#     yield reads, corrupted_reads
#
#
# def corrupt_reads(reads, error_profile, error_loc_rng_rand, base_chose_rng_choice):
#   if len(reads) == 0: return reads
#   reads_per_template = 1 if len(reads[0]) == 1 else 2
#   phred_scores = -10. * numpy.log10(error_profile)
#   qual_str = ''.join([chr(int(min(ep, 93)) + 33) for ep in phred_scores])
#   read_len = phred_scores.size
#
#   # Generate a mirror set of nonsense reads, and then pick from the good read vs bad read depending on a coin toss
#   nonsense_bases = base_chose_rng_choice(['A','C','G','T'], size=len(reads) * reads_per_template * read_len, replace=True, p=[.3, .2, .2, .3])
#   repeated_error_profile = numpy.tile(error_profile, len(reads) * reads_per_template)
#   coin_flip = (error_loc_rng_rand(nonsense_bases.size) < repeated_error_profile).astype('u1')
#
#   corrupted_reads = []
#   for t, template in enumerate(reads):
#     corrupted_template = []
#     for m, read in enumerate(template):
#       idx_st = t * reads_per_template * read_len + read_len * m
#       # Recall a read is a tuple (seq_str, quality_str, coordinate)
#       corrupted_template.append([''.join([base[cf] for cf, base in zip(coin_flip[idx_st:idx_st + read_len], zip(read[0], nonsense_bases[idx_st:idx_st + read_len]))]), qual_str, read[2]])
#     corrupted_reads.append(corrupted_template)
#   return corrupted_reads
#
#
# # 30 us/read
# # def corrupt_reads(reads, error_profile, error_loc_rng_rand, base_chose_rng_choice):
# #   if len(reads) == 0: return reads
# #   reads_per_template = 1 if len(reads[0]) == 1 else 2
# #   phred_scores = -10. * numpy.log10(error_profile)
# #   qual_str = ''.join([chr(int(min(ep, 93)) + 33) for ep in phred_scores])
# #   read_len = phred_scores.size
# #
# #   # Generate a mirror set of nonsense reads, and then pick from the good read vs bad read depending on a coin toss
# #   nonsense_bases = base_chose_rng_choice(['A','C','G','T'], size=len(reads) * reads_per_template * read_len, replace=True, p=[.3, .2, .2, .3])
# #   repeated_error_profile = numpy.tile(error_profile, len(reads) * reads_per_template)
# #   coin_flip = (error_loc_rng_rand(nonsense_bases.size) < repeated_error_profile).astype('u1')
# #
# #   corrupted_reads = []
# #   for t, template in enumerate(reads):
# #     corrupted_template = []
# #     for m, read in enumerate(template):
# #       idx_st = t * reads_per_template * read_len + read_len * m
# #       # Recall a read is a tuple (seq_str, quality_str, coordinate)
# #       corrupted_template.append([''.join([base[cf] for cf, base in zip(coin_flip[idx_st:idx_st + read_len], zip(read[0], nonsense_bases[idx_st:idx_st + read_len]))]), qual_str, read[2]])
# #     corrupted_reads.append(corrupted_template)
# #   return corrupted_reads
#
#
# # # This is inefficently written 46 us/read
# # def corrupt_reads(reads, error_profile, error_loc_rng_rand, base_chose_rng_choice):
# #   if len(reads) == 0: return reads
# #   reads_per_template = 1 if len(reads[0]) == 1 else 2
# #   phred_scores = -10. * numpy.log10(error_profile)
# #   qual_str = ''.join([chr(int(min(ep, 93)) + 33) for ep in phred_scores])
# #   read_len = phred_scores.size
# #   # Generate a mirror set of nonsense reads, and then pick from the good read vs bad read depending on
# #   corrupted_reads = []
# #   for template in reads:
# #     corrupted_template = []
# #     for read in template:
# #       #nonsense_read = base_chose_rng_choice(['A','C','G','T'], size=read_len, replace=True, p=[.3, .2, .2, .3]).tostring()
# #       nonsense_read = base_chose_rng_choice(['A','C','G','T'], size=read_len, replace=True, p=[.3, .2, .2, .3])
# #       coin_flip = (error_loc_rng_rand(read_len) < error_profile).astype('u1')
# #       # Recall a read is a tuple (seq_str, quality_str, coordinate)
# #       corrupted_template.append([''.join([base[cf] for cf, base in zip(coin_flip, zip(read[0], nonsense_read))]), qual_str, read[2]])
# #     corrupted_reads.append(corrupted_template)
# #   return corrupted_reads
#
#
# if __name__ == "__main__":
#   print __explain__
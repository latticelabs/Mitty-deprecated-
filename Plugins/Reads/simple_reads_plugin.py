"""This reads.py plugin generates uniformly sampled reads. It also contains the exponential read corruption model
used by other models (e.g. tiles_reads)

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com

TODO: Fix the corruption algorithm

"""
__explain__ = """
Example parameter file

{
    "model": "simple_reads",
    "args": {
        "paired": false,
        "read_len": 100,
        "template_len": 250,
        "read_loc_rng_seed": 0,
        "read_strand_rng_seed": 1,
        "error_rng_seed": 1,
        "base_chose_rng_seed": 2,
        "max_p_error": 0.8,
        "k": 0.1
    }
}
"""
import numpy
import logging
logger = logging.getLogger(__name__)
# base_sub_mat = {  # GATC For base errors
#                   'G': 'ATC',
#                   'A': 'TCG',
#                   'T': 'CGA',
#                   'C': 'GAT',
# }


def average_read_len(read_len=None, **kwargs):
  """Given the same parameters passed to generate_reads tell us what the average read len is going to be. reads.py
  uses this in combination with coverage and seq_len to figure out how many reads we need."""
  return read_len


def max_read_len(read_len=None, **kwarg):
  """Return maximum read length of reads."""
  return read_len


def read_generator(seq=None,
                   chrom_copy=0,
                   read_start=0,
                   read_stop=0,
                   reads_per_call=1000,
                   num_reads=10000,
                   generate_corrupt_reads=False,
                   paired=False,
                   read_len=None,
                   template_len=None,
                   read_loc_rng_seed=0,
                   read_strand_rng_seed=1,
                   error_rng_seed=1,
                   base_chose_rng_seed=2,
                   max_p_error=.8,
                   k=.1,
                   **kwargs):
  """Given a sequence generate reads with the given characteristics

  Inputs:
    seq              - pair of seq and complement_seq (stringlike) containing the DNA sequence to generate reads from
    chrom_copy       - we use this to change the rng seeds to ensure we take different reads/different positions
                       from each chromosome copy
    read_start       - start generating reads from here (0.0, 1.0)
    read_stop        - stop generating reads from here (0.0, 1.0)
    reads_per_call   - each call to next will generate these many reads
    num_reads        - how many reads do we want in total
    paired           - paired reads or not
    read_len         - Fixed read length (for this model)
    template_len     - Fixed template length (for this model). Only needed if paired is True
    read_loc_rng_seed- Seed for rng that drives the read location picker
    max_p_error      - error probability for last base of read
                       (0.0 is perfect reads, 1.0 -> every base is guaranteed to be wrong)
    k                - exponential factor. 0 < k < 1 The closer this is to 0 the quicker the base error rate drops
    error_loc_rng    - from generate_reads, contains the 2 RNGs we need
    base_chose_rng   - we don't need to return them as they are passed by reference and the state is propagated
                       back up to caller.
    kwargs           - to swallow any other arguments

  Outputs
                                  _________ ( seq_str, quality_str, coordinate)
    perfect_reads  -             /
                 [
                  [[ ... ], [ ... ]],
                  [[ ... ], [ ... ]], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads

  Notes:
  0. This yields a generator
  1. Coordinate is 0-indexed
  2. Quality: Sanger scale 33-126
  3. The number of reads returned on each iteration can be less than reads_per_call because we toss out reads with Ns
     in them

  Read corruption

  1. We use a simple exponential model
  2. We generate random numbers corresponding to base flips all at once to save time.

  """
  assert k <= 1.0
  # We need to initialize the rngs and a bunch of other stuff
  read_loc_rng_randint = numpy.random.RandomState(seed=read_loc_rng_seed + int(chrom_copy)).randint
  read_strand_rng_randint = numpy.random.RandomState(seed=read_strand_rng_seed + int(chrom_copy)).randint
  error_loc_rng_rand = numpy.random.RandomState(seed=error_rng_seed + int(chrom_copy)).rand
  base_chose_rng_choice = numpy.random.RandomState(base_chose_rng_seed + int(chrom_copy)).choice
  error_profile = [max_p_error * k ** n for n in range(read_len)][::-1]

  rl = read_len
  tl = template_len if paired else rl
  if read_start + tl >= read_stop:
    logger.error('Template len too large for given sequence.')
    read_count = num_reads

  if paired:
    num_reads = max(1, num_reads / 2)
    reads_per_call = max(1, reads_per_call / 2)

  read_count = 0
  while read_count < num_reads:
    nominal_read_count = min(reads_per_call, num_reads - read_count)
    rd_st = read_loc_rng_randint(low=read_start, high=read_stop - tl, size=nominal_read_count)
    rd_strand = read_strand_rng_randint(2, size=nominal_read_count)
    reads = []
    if paired:
      for this_rd_st, this_rd_stand in zip(rd_st, rd_strand):
        seq_1 = seq[this_rd_stand][this_rd_st:this_rd_st + rl]
        seq_2 = seq[1 - this_rd_stand][this_rd_st + tl - 1:this_rd_st + tl - 1 - rl:-1]
        if 'N' in seq_1 or 'N' in seq_2:  # read taken from a masked/unknown region
          continue
        reads.append([[seq_1, '~' * rl, this_rd_st],
                      [seq_2, '~' * rl, this_rd_st + tl - rl]])
    else:
      for this_rd_st in rd_st:
        if 'N' in seq[this_rd_st:this_rd_st + tl]:  # read taken from a masked/unknown region
          continue
        reads.append([[seq[this_rd_st:this_rd_st + rl], '~' * rl, this_rd_st]])

    if generate_corrupt_reads:
      corrupted_reads = corrupt_reads(reads, error_profile, error_loc_rng_rand, base_chose_rng_choice)
    else:
      corrupted_reads = None

    read_count += nominal_read_count * (2 if paired else 1) # len(reads) We don't use actual read count to avoid thrashing when we have tons of Ns in the sequence
    logger.debug('Generated {:d} reads'.format(read_count))
    yield reads, corrupted_reads


def corrupt_reads(reads, error_profile, error_loc_rng_rand, base_chose_rng_choice):
  if len(reads) == 0: return reads
  phred_scores = [-10. * numpy.log10(ep) for ep in error_profile]
  qual_str = ''.join([chr(int(min(ep, 93)) + 33) for ep in phred_scores])
  read_len = len(reads[0][0][0])

  # Generate a mirror set of nonsense reads, and then pick from the good read vs bad read depending on
  corrupted_reads = []
  for template in reads:
    corrupted_template = []
    for read in template:
      nonsense_read = base_chose_rng_choice(['A','C','G','T'], size=read_len, replace=True, p=[.3, .2, .2, .3]).tostring()
      coin_flip = (error_loc_rng_rand(read_len) < error_profile).astype('u1')
      # Recall a read is a tuple (seq_str, quality_str, coordinate)
      corrupted_template.append([''.join([base[cf] for cf, base in zip(coin_flip, zip(read[0], nonsense_read))]), qual_str, read[2]])
    corrupted_reads.append(corrupted_template)
  return corrupted_reads



# def corrupt_reads(reads, read_len=100,
#                   error_rng_seed=1,
#                   base_chose_rng_seed=2,
#                   max_p_error=.8,
#                   k=.1,
#                   **kwargs):
#   """Simple exponential error model for reads.
#   Inputs:
#     reads        - perfect reads as produced by the read_generator
#     read_len     - the (fixed) length of the reads
#     max_p_error  - error probability for last base of read
#                    (0.0 is perfect reads, 1.0 -> every base is guaranteed to be wrong)
#     k            - exponential factor. 0 < k < 1 The closer this is to 0 the quicker the base error rate drops
#     error_loc_rng  - from generate_reads, contains the 2 RNGs we need
#     base_chose_rng - we don't need to return them as they are passed by reference and the state is propagated
#                      back up to caller.
#
#   Outputs:
#                                   _________ ( seq_str, quality_str, coordinate)
#     corrupt_reads  -             /
#                  [
#                   [[ ... ], [ ... ]],
#                   [[ ... ], [ ... ]], -> inner list = 2 elements if paired reads, 1 otherwise
#                        .
#                        .
#                        .
#                  ] -> outer list = number of reads
#
#   Notes:
#   1. The only tricky thing about this function is that it has a state (stored as an attribute): We need to carry the
#   state of the random number generators from call to call. We use try: except: for the initialization of the generators
#   as this is slightly faster than if hasattr ... You will note that read_generator also stores the state of its rngs,
#   but because it is structured as a generator, we don't think twice about this.
#
#   """
#   if len(reads) == 0: return reads  # We return `reads` so we match the shape (paired or not paired)
#
#
#   # A little sleight of hand here. We store state in a function attribute. The only caveat here is that we can't use
#   # reuse this function expecting to be able to store separate states.
#   # It is faster to ask for forgiveness, than permission
#   # (http://assorted-experience.blogspot.com/2014/05/storing-state-in-python-function.html)
#   try:
#     base_errors = corrupt_reads.error_loc_rng.rand(len(reads[0]), len(reads), read_len)  # Coin toss to see if we error the base call
#     base_subs = corrupt_reads.base_chose_rng.randint(3, size=(len(reads[0]), len(reads), read_len))  # If so, what base will be call it
#     error_profile = corrupt_reads.error_profile
#   except AttributeError:
#     corrupt_reads.error_loc_rng = numpy.random.RandomState(seed=error_rng_seed)
#     corrupt_reads.base_chose_rng = numpy.random.RandomState(base_chose_rng_seed)
#     corrupt_reads.error_profile = [max_p_error * k ** n for n in range(read_len)][::-1]
#
#     base_errors = corrupt_reads.error_loc_rng.rand(len(reads[0]), len(reads), read_len)  # Coin toss to see if we error the base call
#     base_subs = corrupt_reads.base_chose_rng.randint(3, size=(len(reads[0]), len(reads), read_len))  # If so, what base will be call it
#
#   rev_error_profile = [max_p_error * k ** n for n in range(read_len)]  # For the second of the pair
#
#
#
#
#   rev_error_profile = [max_p_error * k ** n for n in range(read_len)]  # For the second of the pair
#   error_profile = [rev_error_profile[::-1], rev_error_profile]
#   # error_profile[0] -> base call error curve for forward read (5' end more error prone)
#   # error_profile[1] -> base call error curve for backward read (3' end more error prone)
#   qual = [''.join([chr(int(126-(126-33)*ep)) for ep in error_profile[n]]) for n in [0,1]]
#   corrupted_reads = []
#
#   for n in range(len(reads)):
#     these_reads = [[None]]*len(reads[0])
#     for p in range(len(reads[0])):
#       this_read = bytearray(reads[n][p][0])
#       for m in range(read_len):
#         if base_errors[0, n, m] < error_profile[p][m]:
#           ref = reads[n][p][0][m]
#           this_read[m] = base_sub_mat.get(ref,'NNN')[base_subs[0, n, m]]  # The 'NNN' ensures we return an N for
#                                                                           # any ambiguous bases
#       these_reads[p] = [this_read.__str__(), qual[p], reads[n][p][2]]
#     corrupted_reads.append(these_reads)
#
#   return corrupted_reads

if __name__ == "__main__":
  import sys
  if len(sys.argv) == 2:  # Print explain
    if sys.argv[1] == 'explain':
      print __explain__
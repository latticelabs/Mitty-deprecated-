"""This reads.py plugin generates uniformly sampled reads. It also contains the exponential read corruption model
used by other models (e.g. tiles_reads)

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
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
base_sub_mat = {  # GATC
                  'G': 'ATC',
                  'A': 'TCG',
                  'T': 'CGA',
                  'C': 'GAT',
}


def average_read_len(read_len=None, **kwargs):
  """Given the same parameters passed to generate_reads tell us what the average read len is going to be. reads.py
  uses this in combination with coverage and seq_len to figure out how many reads we need."""
  return read_len


def read_generator(seq,
                   read_start=0,
                   read_stop=0,
                   reads_per_call=1000,
                   num_reads=10000,
                   paired=False,
                   read_len=None,
                   template_len=None,
                   read_loc_rng_seed=0,
                   **kwargs):
  """Given a sequence generate reads with the given characteristics

  Inputs:
    seq              - string(like)s containing the DNA sequence to generate reads from
    seq_len          - length of the sequence
    read_start       - start generating reads from here (0.0, 1.0)
    read_stop        - stop generating reads from here (0.0, 1.0)
    reads_per_call   - each call to next will generate these many reads
    num_reads        - how many reads do we want in total
    paired           - paired reads or not
    read_len         - Fixed read length (for this model)
    template_len     - Fixed template length (for this model). Only needed if paired is True
    read_loc_rng_seed- Seed for rng that drives the read location picker
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
  """
  # We need to initialize the rngs and a bunch of other stuff
  read_loc_rng = numpy.random.RandomState(seed=read_loc_rng_seed)

  rl = read_len
  tl = template_len if paired else rl
  read_count = 0
  if read_start + tl >= read_stop:
    logger.error('Template len too large for given sequence.')
    read_count = num_reads

  if paired:
    num_reads = max(1, num_reads / 2)
    reads_per_call = max(1, reads_per_call / 2)

  while read_count < num_reads:
    nominal_read_count = min(reads_per_call, num_reads - read_count)
    rd_st = read_loc_rng.randint(low=read_start, high=read_stop - tl, size=nominal_read_count)
    reads = []
    if paired:
      for this_rd_st in rd_st:
        if 'N' in seq[this_rd_st:this_rd_st + tl]:  # read taken from a masked/unknown region
          continue
        reads.append([[seq[this_rd_st:this_rd_st + rl], '~' * rl, this_rd_st],
                      [seq[this_rd_st + tl - rl:this_rd_st + tl], '~' * rl, this_rd_st + tl - rl]])
    else:
      for this_rd_st in rd_st:
        if 'N' in seq[this_rd_st:this_rd_st + tl]:  # read taken from a masked/unknown region
          continue
        reads.append([[seq[this_rd_st:this_rd_st + rl], '~' * rl, this_rd_st]])

    read_count += nominal_read_count * (2 if paired else 1) # len(reads) We don't use actual read count to avoid thrashing when we have tons of Ns in the sequence
    yield reads


def corrupt_reads(reads, read_len=100,
                  error_rng_seed=1,
                  base_chose_rng_seed=2,
                  max_p_error=.8,
                  k=.1,
                  **kwargs):
  """Simple exponential error model for reads.
  Inputs:
    reads        - perfect reads as produced by the read_generator
    read_len     - the (fixed) length of the reads
    max_p_error  - error probability for last base of read
                   (0.0 is perfect reads, 1.0 -> every base is guaranteed to be wrong)
    k            - exponential factor. 0 < k < 1 The closer this is to 0 the quicker the base error rate drops
    error_loc_rng  - from generate_reads, contains the 2 RNGs we need
    base_chose_rng - we don't need to return them as they are passed by reference and the state is propagated
                     back up to caller.

  Outputs:
                                  _________ ( seq_str, quality_str, coordinate)
    corrupt_reads  -             /
                 [
                  [[ ... ], [ ... ]],
                  [[ ... ], [ ... ]], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads

  Note:
  The only tricky thing about this function is that it has a state (stored as an attribute): We need to carry the
  state of the random number generators from call to call. We use try: except: for the initialization of the generators
  as this is slightly faster than if hasattr ... You will note that read_generator also stores the state of its rngs,
  but because it is structured as a generator, we don't think twice about this.
  """
  if len(reads) == 0: return reads  # We return `reads` so we match the shape (paired or not paired)
  rev_error_profile = [max_p_error * k ** n for n in range(read_len)]  # For the second of the pair
  error_profile = [rev_error_profile[::-1], rev_error_profile]
  # error_profile[0] -> base call error curve for forward read (5' end more error prone)
  # error_profile[1] -> base call error curve for backward read (3' end more error prone)
  qual = [''.join([chr(int(126-(126-33)*ep)) for ep in error_profile[n]]) for n in [0,1]]
  corrupted_reads = []
  # A little sleight of hand here. We store state in a function attribute. The only caveat here is that we can't use
  # reuse this function expecting to be able to store separate states.
  try:
    base_errors = corrupt_reads.error_loc_rng.rand(len(reads[0]), len(reads), read_len)  # Coin toss to see if we error the base call
    base_subs = corrupt_reads.base_chose_rng.randint(3, size=(len(reads[0]), len(reads), read_len))  # If so, what base will be call it
  # It is faster to ask for forgiveness, than permission
  # (http://assorted-experience.blogspot.com/2014/05/storing-state-in-python-function.html)
  except AttributeError:
    corrupt_reads.error_loc_rng = numpy.random.RandomState(seed=error_rng_seed)
    corrupt_reads.base_chose_rng = numpy.random.RandomState(base_chose_rng_seed)
    base_errors = corrupt_reads.error_loc_rng.rand(len(reads[0]), len(reads), read_len)  # Coin toss to see if we error the base call
    base_subs = corrupt_reads.base_chose_rng.randint(3, size=(len(reads[0]), len(reads), read_len))  # If so, what base will be call it

  for n in range(len(reads)):
    these_reads = [[None]]*len(reads[0])
    for p in range(len(reads[0])):
      this_read = bytearray(reads[n][p][0])
      for m in range(read_len):
        if base_errors[0, n, m] < error_profile[p][m]:
          ref = reads[n][p][0][m]
          this_read[m] = base_sub_mat.get(ref,'NNN')[base_subs[0, n, m]]  # The 'NNN' ensures we return an N for
                                                                          # any ambiguous bases
      these_reads[p] = [this_read.__str__(), qual[p], reads[n][p][2]]
    corrupted_reads.append(these_reads)

  return corrupted_reads

if __name__ == "__main__":
  import sys
  if len(sys.argv) == 2:  # Print explain
    if sys.argv[1] == 'explain':
      print __explain__
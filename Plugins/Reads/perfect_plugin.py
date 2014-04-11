import numpy
import logging
logger = logging.getLogger(__name__)
base_sub_mat = {  # GATC
                  'G': 'ATC',
                  'A': 'TCG',
                  'T': 'CGA',
                  'C': 'GAT'
}


def average_read_len(read_len, **kwargs):
  """Given the same parameters passed to generate_reads tell us what the average read len is going to be"""
  return read_len


def generate_reads(seq,
                   start_reads=0,
                   stop_reads=0,
                   num_reads=1000,
                   paired=False,
                   read_len=None,
                   template_len=None,
                   read_loc_rng_seed=0,
                   error_rng_seed=1,
                   base_chose_rng_seed=2,
                   max_p_error=.8,
                   k=.1,
                   prev_state=None,
                   **kwargs):
  """Given a list of sequences generate reads with the given characteristics

  Inputs:
    seq              - string(like)s containing the DNA sequence to generate reads from
    start            - 0-indexed coordinate of start of section reads will be generated from
    stop             - 0-indexed coordinate of end of section reads will be generated from (-1 means till end)
    num_reads        - reads to generate this call to the function
    paired           - paired reads or not
    read_len         - Fixed read length
    template_len     - Template length. Only needed if paired is True
    read_loc_rng_seed- Seed for rng that drives the read location picker
    error_rng_seed   - Seed for rng that determines if there is a read error or not on the base
    base_chose_rng_seed - Seed for rng that determines which base is erroneously read
    prev_state       - previous state carried over from call to call.
                       In this case it stores the three rngs
    kwargs           - to swallow any other arguments

  Outputs
                                 _________ ( seq_str, quality_str, coordinate)
    reads     -  [              /
                  [( ... ), ( ...)],
                  [( ... ), ( ...)], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads

    perfect_reads - same format as reads, same size, but with no read errors

  Quality: Sanger scale 33-126
  """
  if prev_state is None:  # If we are running this for the first time, we need to initialize the rngs
    read_loc_rng = numpy.random.RandomState(seed=read_loc_rng_seed)
    error_loc_rng = numpy.random.RandomState(seed=error_rng_seed)
    base_chose_rng = numpy.random.RandomState(base_chose_rng_seed)
  else:
    read_loc_rng = prev_state['read_loc_rng']
    error_loc_rng = prev_state['error_loc_rng']
    base_chose_rng = prev_state['base_chose_rng']

  rl = read_len
  tl = template_len
  if (stop_reads - tl < read_len) and paired:
    logger.error('Template len too large for given sequence.')

  if paired:
    rd_st = read_loc_rng.randint(low=start_reads, high=stop_reads - tl, size=num_reads)  # Reads are uniformly distributed
    reads = [[(seq[rd_st[n]:rd_st[n] + rl], '~' * rl, rd_st[n]+1),
              (seq[rd_st[n] + tl - rl:rd_st[n] + tl], '~' * rl, rd_st[n] + tl - rl+1)] for n in range(num_reads)]
  else:
    rd_st = read_loc_rng.randint(low=start_reads, high=stop_reads - rl, size=num_reads)
    reads = [[(seq[rd_st[n]:rd_st[n] + rl], '~' * rl, rd_st[n]+1)] for n in range(num_reads)]

  corr_reads = corrupt_reads_expon(reads, read_len, max_p_error, k, error_loc_rng, base_chose_rng)

  prev_state = {
    'read_loc_rng': read_loc_rng,
    'error_loc_rng': error_loc_rng,
    'base_chose_rng': base_chose_rng
  }

  return corr_reads, reads, prev_state


def corrupt_reads_expon(reads, read_len=100, max_p_error=.8, k=.1, error_loc_rng=None, base_chose_rng=None):
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
    reads        - corrupted reads
  """
  if len(reads) == 0: return reads  # We return `reads` so we match the shape (paired or not paired)
  rev_error_profile = [max_p_error * k ** n for n in range(read_len)]  # For the second of the pair
  error_profile = [rev_error_profile[::-1], rev_error_profile]
  # error_profile[0] -> base call error curve for forward read (5' end more error prone)
  # error_profile[1] -> base call error curve for backward read (3' end more error prone)
  qual = [''.join([chr(int(126-(126-33)*ep)) for ep in error_profile[n]]) for n in [0,1]]
  corrupted_reads = []
  base_errors = error_loc_rng.rand(len(reads[0]), len(reads), read_len)  # Coin toss to see if we error the base call
  base_subs = base_chose_rng.randint(3, size=(len(reads[0]), len(reads), read_len))  # If so, what base will be call it
  for n in range(len(reads)):
    these_reads = [[None]]*len(reads[0])
    for p in range(len(reads[0])):
      this_read = bytearray(reads[n][p][0])
      for m in range(read_len):
        if base_errors[0, n, m] < error_profile[p][m]:
          ref = reads[n][p][0][m]
          this_read[m] = base_sub_mat[ref][base_subs[0, n, m]]
      these_reads[p] = [this_read.__str__(), qual[p], reads[n][p][2]]
    corrupted_reads.append(these_reads)

  return corrupted_reads
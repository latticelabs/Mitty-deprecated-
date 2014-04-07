import numpy
import logging

logger = logging.getLogger(__name__)


def generate_reads(seq,
                   start=0,
                   stop=-1,
                   num_reads=1000,
                   read_len=None,
                   template_len=None,
                   paired=False,
                   rng_seed=0,
                   prev_state=None):
  """Given a list of sequences generate reads with the given characteristics

  Inputs:
    seq              - string(like)s containing the DNA sequence to generate reads from
    num_reads        - reads to generate this call to the function
    read_len         - Fixed read length
    template_len     - Template length. Only needed if paired is True
    coverage         - the coverage that we want
    paired           - paired reads or not
    rng              - numpy.random.RandomState(seed)

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
  if prev_state is None:  # If we are running this for the first time, we need to initialize the rng
    rng = numpy.random.RandomState(seed=rng_seed)
  else:
    rng = prev_state

  if stop==-1: stop=len(seq)
  logger.debug('Starting to generate reads')
  rl = read_len
  tl = template_len
  if (stop - tl < read_len) and paired:
    logger.error('Template len too large for given sequence.')

  if paired:
    rd_st = rng.randint(low=start, high=stop - tl, size=num_reads)  # Reads are uniformly distributed
    reads = [[(seq[rd_st[n]:rd_st[n] + rl], '~' * rl, rd_st[n]),
              (seq[rd_st[n] + tl - rl:rd_st[n] + tl], '~' * rl, rd_st[n] + tl - rl)] for n in range(num_reads)]
  else:
    rd_st = rng.randint(low=start, high=stop - rl, size=num_reads)
    reads = [[(seq[rd_st[n]:rd_st[n] + rl], '~' * rl, rd_st[n])] for n in range(num_reads)]
  logger.debug('Finished generating reads')
  return reads, reads, rng

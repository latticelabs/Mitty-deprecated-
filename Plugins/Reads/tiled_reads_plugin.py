"""This reads.py plugin generates tiled reads with a guaranteed coverage for every base. The reads are taken
deterministically, starting at the beginning of the sequence and then taking reads with the requested overlap.
At the end of the sequence the overlap is adjusted to make sure every base is covered.

Example parameter file
{
    "__comment__": "Example parameter file for reads program. Seven Bridges Genomics Current contact: kaushik.ghose@sbgenomics.com",
    "model_name": "tiled_reads",
    "seq_header": "gi|52547303|gb|AY735451.1| Porcine circovirus isolate Hebei capsid protein gene, complete cds",
    "seq_file": "Data/porcine_circovirus.smalla",
    "corrupted_reads_file": "Data/corrupted_reads.bam",
    "perfect_reads_file": "Data/perfect_reads.bam",
    "output_type": "bam",
    "start": 0,
    "stop": -1,
    "coverage": 2,
    "average_read_len": 100,
    "args": {
        "paired": true,
        "read_len": 100,
        "template_len": 250,
        "overlap": 0.5
    }
}

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
import logging

logger = logging.getLogger(__name__)


def generate_reads(seq,
                   start=0,
                   stop=-1,
                   num_reads=1000,
                   paired=False,
                   read_len=100,
                   template_len=250,
                   overlap=0.5,
                   prev_state=None):
  """Given a list of sequences generate reads with the given characteristics

  Inputs:
    seq              - string(like)s containing the DNA sequence to generate reads from
    start            - 0-indexed coordinate of start of section reads will be generated from
    stop             - 0-indexed coordinate of end of section reads will be generated from (-1 means till end)
    num_reads        - reads to generate this call to the function
    paired           - paired reads or not
    read_len         - Fixed read length
    template_len     - Template length. Only needed if paired is True
    overlap          - what fraction of the previous read should this read cover
    prev_state       - previous state carried over from call to call.
                       In this case it is simply start position of the next read we should generate

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
    next_read_start = start
  else:
    next_read_start = prev_state
  if stop==-1: stop=len(seq)
  if (stop - template_len < read_len) and paired:
    logger.error('Template len too large for given sequence.')

  logger.debug('Starting to generate reads')
  nrs = next_read_start
  rl = read_len
  tl = template_len if paired else rl
  ra = int((1. - overlap) * read_len)
  # TODO: Write this more Pythonically?
  rd_st = nrs
  reads = []
  for n in range(num_reads):
    rd_e = rd_st + tl
    if rd_e > stop:  # We need to slide back to make sure we cover the end of the sequence
      rd_e = stop
      rd_st = rd_e - tl
    this_read = [(seq[rd_st:rd_st+rl], '~' * rl, rd_st + 1)] # Coordinates are 1 indexed so we need +1
    if paired:
      this_read += [(seq[rd_e-rl:rd_e], '~' * rl, rd_e - rl + 1)]
    reads += [this_read]

    if rd_e == stop:  # We just covered the sequence once, we need to reset
      rd_st = 0
    else:  # Advance as usual
      rd_st += ra

  logger.debug('Finished generating reads')
  return reads, reads, rd_st

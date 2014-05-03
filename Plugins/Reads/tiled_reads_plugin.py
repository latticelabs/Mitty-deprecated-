"""This reads.py plugin generates tiled reads with a guaranteed coverage for every base. The reads are taken
deterministically, starting at the beginning of the sequence and then taking reads with the requested overlap. The
overlap of the last read is adjusted so we cover every base. When we get to the end of the sequence we wrap around and
start with a shifted offset.

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
import numpy
import logging
logger = logging.getLogger(__name__)

def average_read_len(read_len, **kwargs):
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
                   read_advance=None,
                   **kwargs):
  """Given a sequence generate reads with the given characteristics

  Inputs:
    seq              - string(like)s containing the DNA sequence to generate reads from
    read_start       - start generating reads from here
    read_stop        - stop generating reads from here
    reads_per_call   - each call to next will generate these many reads
    num_reads        - how many reads do we want in total
    paired           - paired reads or not
    read_len         - Fixed read length (for this model)
    template_len     - Fixed template length (for this model). Only needed if paired is True
    read_advance     - Advance these many bases to take the next read
    kwargs           - to swallow any other arguments

  Outputs
    perfect_reads    -

                                 _________ ( seq_str, quality_str, coordinate)
                                /
                 [
                  [( ... ), ( ...)],
                  [( ... ), ( ...)], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads

  Note: coordinate is 0-indexed

  Quality: Sanger scale 33-126
  """
  rl = read_len
  tl = template_len if paired else rl
  read_count = 0
  if read_start + tl >= read_stop:
    logger.error('Template len too large for given sequence.')
    read_count = num_reads

  read_offset = 0
  nominal_read_start = read_start
  while read_count < num_reads:
    reads = []
    for n in range(min(reads_per_call, num_reads-read_count)):
      if nominal_read_start + tl > read_stop:
        this_read_start = read_stop - tl
        read_offset += 1
        if read_offset > rl:
          read_offset = 0
        nominal_read_start = read_start + read_offset
      else:
        this_read_start = nominal_read_start
        nominal_read_start += read_advance

      these_reads = [[seq[this_read_start:this_read_start + rl], '~' * rl, this_read_start]]
      if paired:
        these_reads += [[seq[this_read_start + tl - rl:this_read_start + tl], '~' * rl, this_read_start + tl - rl]]
      reads.append(these_reads)

    read_count += len(reads)
    yield reads
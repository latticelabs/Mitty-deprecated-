"""A dummy base class for read plugins to simplify some things"""
import numpy as np


class ReadModel:
  """Base class for read plugins"""
  dtype = [('start_a', 'i4'), ('read_len', 'i4'), ('read_order', 'i1'),
           ('perfect_reads', 'O'), ('corrupt_reads', 'O'), ('phred', 'O')]

  def __init__(self, paired):
    self.paired = paired

  def get_zero_reads(self):
    """Return empty array of reads. Useful for concatenation etc."""
    return np.recarray(dtype=ReadModel.dtype, shape=0), self.paired

  def get_reads(self, seq, seq_c, start_base=0, end_base=None, coverage=0.01, corrupt=False, seed=1):
    """The main simulation calls this function.

    :param seq:      forward sequence
    :param seq_c:    complement of sequence
    :param start_base: base to start taking reads from
    :param end_base: base to stop
    :param coverage: coverage
    :param corrupt:  T/F whether we should compute corrupted read or not
    :param seed:     random number generator seed
    :return: reads, paired

    reads is a numpy recarray with the following fields
      'start_a'  -> start index on the seq
      'read_len' -> length of the read
      'read_order' ->  0 or 1, indicating if the read is from the forward or reverse strand
      'perfect_reads'     -> perfect read sequence
      'corrupt_reads'   -> corrupted read sequence
      'phred'  -> phred score string for base quality

    paired indicates if the reads are in pairs or not
    """
    return self.get_zero_reads()
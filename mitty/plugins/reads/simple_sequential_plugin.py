"""Very simple deterministic read generator with no read corruption.

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
__example_param_text = """
{
  'read_len': 100,          # length of each read
  'template_len': 250,      # length of template, only needed for paired reads
  'paired': True            # Give paired reads?
}
"""

_description = """
This is a simple, deterministic read generator with no read corruption.
Example parameter set:
""" + __example_param_text

_example_params = eval(__example_param_text)

import numpy as np

from mitty.plugins.reads.base_plugin import ReadModel

import logging
logger = logging.getLogger(__name__)


class Model(ReadModel):
  def __init__(self, read_len=100, template_len=250, paired=True):
    """."""
    self.read_len, self.template_len = read_len, template_len
    self.start_base = 0
    ReadModel.__init__(self, paired)

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
      'rev_complement' ->  0 or 1, indicating if the read is from the forward or reverse strand
      'perfect_read'     -> perfect read sequence
      'corrupted_read'   -> corrupted read sequence

    paired indicates if the reads are in pairs or not
    """
    assert len(seq) > self.template_len, 'Template size should be less than sequence length'
    while seq[self.start_base] == 'N' or self.start_base < start_base:
      self.start_base += 10

    stride = float(self.read_len) / coverage
    template_locs = np.array([x for x in np.arange(self.start_base, end_base or len(seq), stride, dtype='i4') if seq[x] != 'N'], dtype='i4')

    reads = np.core.recarray(dtype=ReadModel.dtype, shape=(2 if self.paired else 1) * template_locs.shape[0])
    reads['read_len'] = self.read_len
    if self.paired:
      reads['start_a'][::2] = template_locs
      reads['start_a'][1::2] = template_locs + self.template_len - self.read_len
      reads['read_order'][::2] = 0
      reads['read_order'][1::2] = 1
    else:
      reads['start_a'] = template_locs
      reads['read_len'] = self.read_len
      reads['read_order'] = 0

    r_start, r_len, r_o, pr, cr = reads['start_a'], reads['read_len'], reads['read_order'], reads['perfect_reads'], reads['corrupt_reads']
    for n in xrange(reads.shape[0]):
      if r_o[n] == 0:  # Forward read
        pr[n] = seq[r_start[n]:r_start[n] + r_len[n]]
      else:  # Reverse strand read
        pr[n] = seq_c[r_start[n]:r_start[n] + r_len[n]][::-1]
      if corrupt:
        cr[n] = pr[n]

    return reads, self.paired


def self_test():
  """Basic self test"""
  import string
  DNA_complement = string.maketrans('ATCGN', 'TAGCN')
  seq = 'ATGTCGCCGGGCGCCATGCGTGCCGTTGTTCCCATTATCCCATTCCTTTTGGTTCTTGTCGGTGTATCGGGGGTTCCCACCAACGTCTCCTCCACCACCCAACCCCAACTCCAGACCACCGGTCGTCCCTCGCATGAAGCCCCCAACATGACCCAGACCGGCACCACCGACTCTCCCACCGCCATCAGCCTTACCACGCCCGACCACACACCCCCCATGCCAAGTATCGGACTGGAGGAGGAGGAAGAGGAGGAGGGGGCCGGGGATGGCGAACATCTTGAGGGGGGAGATGGGACCCGTGACACCCTACCCCAGTCCCCGGGTCCAGCCGTCCCGTTGGCCGGGGATGACGAGAAGGACAAACCCAACCGTCCCGTAGTCCCACCCCCCGGTCCCAACAACTCCCCCGCGCGCCCCGAGACCAGTCGACCGAAGACACCCCCCACCAGTATCGGGCCGCTGGCAACTCGACCCACGACCCAACTCCCCTCAAAGGGGCGACCCTTGGTTCCGACGCCTCAACATACCCCGCTGTTCTCGTTCCTCACTGCCTCCCCCGCCCTGGACACCCTCTTCGTCGTCAGCACCGTCATCCACACCTTATCGTTTTTGTGTATTGTTGCGATGGCGACACACCTGTGTGGCGGTTGGTCCAGACGCGGGCGACGCACACACCCTAGCGTGCGTTACGTGTGCCTGCCGCCCGAACGCGGGTAG'
  seq_c = string.translate(seq, DNA_complement)
  mdl = Model(4, 8, True)
  rd, paired = mdl.get_reads(seq, seq_c, start_base=0, end_base=len(seq), coverage=.00001, corrupt=True)
  assert type(rd) == np.core.records.recarray  # Basically, the previous code should just run
  assert paired is True


if __name__ == "__main__":
  print _description
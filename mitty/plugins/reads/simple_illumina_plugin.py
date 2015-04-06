"""This is the stock read plugin that approximates illumina reads

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
__example_param_text = """
{
  'read_len': 100,          # length of each read
  'template_len_mean': 250, # mean length of template
  'template_len_sd': 30,    # sd of template length
  'max_p_error': 0.01,      # Maximum error rate at tip of read
  'k': 20                   # exponent
}
"""

_description = """
This read generator generates Illumina like reads with a exponential error profile
Example parameter set:
""" + __example_param_text

_example_params = eval(__example_param_text)

import numpy as np

import mitty.lib.util as mutil

import logging
logger = logging.getLogger(__name__)


class Model:
  def __init__(self, read_len=100, template_len_mean=250, template_len_sd=50, max_p_error=0.01, k=20):
    """."""
    self.read_len, self.template_len_mean, self.template_len_sd = read_len, template_len_mean, template_len_sd
    l = np.linspace(1e-10, 1, self.read_len)
    self.error_profile = max_p_error * (np.exp(k*l) - 1)/(np.exp(k) - 1)
    self.phred = ''.join([chr(int(33 + max(0, min(-10*np.log10(p), 93)))) for p in self.error_profile])

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
    assert len(seq) > self.template_len_mean * 3, 'Template size should be less than 1/3rd sequence length'
    p_template = 0.5 * coverage / float(self.read_len)  # Per base probability of a read
    template_loc_rng, read_order_rng, template_len_rng, error_loc_rng, base_choice_rng = mutil.initialize_rngs(seed, 5)

    template_locs = np.array([x for x in mutil.place_poisson(template_loc_rng, p_template, start_base, min(end_base or len(seq), len(seq) - self.template_len_mean * 3)) if seq[x] != 'N'], dtype='i4')
    # TODO, refactor this into mutil by having seq as an numpy array rather than string
    # also change for variant plugins. Need to change io.py. Should be transparent (non breaking change) at io.py
    template_lens = (template_len_rng.randn(template_locs.shape[0]) * self.template_len_sd + self.template_len_mean).astype('i4')
    np.clip(template_lens, self.read_len, self.template_len_mean * 3, template_lens)
    read_order = read_order_rng.randint(2, size=template_locs.shape[0])  # Which read comes first?

    dtype = [('start_a', 'i4'), ('read_len', 'i4'), ('read_order', 'i1'),
             ('perfect_reads', 'S' + str(self.read_len)), ('corrupt_reads', 'S' + str(self.read_len)),
             ('phred', 'S' + str(self.read_len))]
    reads = np.core.recarray(dtype=dtype, shape=2 * template_locs.shape[0])
    reads['start_a'][::2] = template_locs
    reads['start_a'][1::2] = template_locs + template_lens - self.read_len
    reads['read_len'] = self.read_len
    reads['read_order'][::2] = read_order[:]
    reads['read_order'][1::2] = 1 - read_order[:]
    r_start, r_len, r_o, pr = reads['start_a'], reads['read_len'], reads['read_order'], reads['perfect_reads']
    for n in xrange(reads.shape[0]):
      if r_o[n] == 0:  # Forward read
        pr[n] = seq[r_start[n]:r_start[n] + r_len[n]]
      else:  # Reverse strand read
        pr[n] = seq_c[r_start[n]:r_start[n] + r_len[n]][::-1]
    if corrupt:
      self.corrupt_reads(reads, error_loc_rng, base_choice_rng)
    return reads, True

  def corrupt_reads(self, reads, error_loc_rng, base_choice_rng):
    """Corrupt reads

    :param reads:   with the start_a, read_len and read_order fields filled out
    :return: Fill in corrupted reads in place
    """
    pr, cr, ph = reads['perfect_reads'], reads['corrupt_reads'], reads['phred']
    p = error_loc_rng.rand(reads.shape[0], self.read_len)
    idx = [np.where(p[n, :] < self.error_profile)[0] for n in xrange(reads.shape[0])]
    idx_tot = 0
    for n in xrange(reads.shape[0]):
      idx_tot += idx[n].size
    corrupt_bases = base_choice_rng.choice(['A','C','G','T'], size=idx_tot, replace=True, p=[.3, .2, .2, .3]).tostring()

    offset = 0
    for n in xrange(reads.shape[0]):
      ph[n] = self.phred
      #idx, = np.where(p[n, :] < self.error_profile)
      #corrupt_bases = base_choice_rng.choice(['A','C','G','T'], size=idx.size, replace=True, p=[.3, .2, .2, .3]).tostring()
      b = bytearray(pr[n])
      for m, i in enumerate(idx[n]):
        b[i] = corrupt_bases[offset +m]
      cr[n] = str(b)
      offset += idx[n].size


def self_test():
  """Basic self test"""
  import string
  DNA_complement = string.maketrans('ATCGN', 'TAGCN')
  seq = 'ATGTCGCCGGGCGCCATGCGTGCCGTTGTTCCCATTATCCCATTCCTTTTGGTTCTTGTCGGTGTATCGGGGGTTCCCACCAACGTCTCCTCCACCACCCAACCCCAACTCCAGACCACCGGTCGTCCCTCGCATGAAGCCCCCAACATGACCCAGACCGGCACCACCGACTCTCCCACCGCCATCAGCCTTACCACGCCCGACCACACACCCCCCATGCCAAGTATCGGACTGGAGGAGGAGGAAGAGGAGGAGGGGGCCGGGGATGGCGAACATCTTGAGGGGGGAGATGGGACCCGTGACACCCTACCCCAGTCCCCGGGTCCAGCCGTCCCGTTGGCCGGGGATGACGAGAAGGACAAACCCAACCGTCCCGTAGTCCCACCCCCCGGTCCCAACAACTCCCCCGCGCGCCCCGAGACCAGTCGACCGAAGACACCCCCCACCAGTATCGGGCCGCTGGCAACTCGACCCACGACCCAACTCCCCTCAAAGGGGCGACCCTTGGTTCCGACGCCTCAACATACCCCGCTGTTCTCGTTCCTCACTGCCTCCCCCGCCCTGGACACCCTCTTCGTCGTCAGCACCGTCATCCACACCTTATCGTTTTTGTGTATTGTTGCGATGGCGACACACCTGTGTGGCGGTTGGTCCAGACGCGGGCGACGCACACACCCTAGCGTGCGTTACGTGTGCCTGCCGCCCGAACGCGGGTAG'
  seq_c = string.translate(seq, DNA_complement)
  mdl = Model(4, 8, 2, max_p_error=1)
  rd, paired = mdl.get_reads(seq, seq_c, start_base=0, end_base=len(seq), coverage=.00001, corrupt=True)
  assert type(rd) == np.core.records.recarray  # Basically, the previous code should just run
  assert paired == True

if __name__ == "__main__":
  print _description
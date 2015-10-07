"""This is the stock read plugin that approximates Illumina reads

GC bias model based on:

"Summarizing and correcting the GC content bias in high-throughput sequencing", Benjamini and Speed,
Nucleic Acids Research, 2012
http://nar.oxfordjournals.org/content/early/2012/02/08/nar.gks001.full

1. "GC effect is mostly driven by the GC composition of the full fragment."
2. "the GC curve is unimodal is key to this analysis. In all data sets shown,
the rate of GC-poor or GC-rich fragments is significantly lower than average, in many cases zero."

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
__example_param_text = """
{
  'read_len': 100,          # length of each read
  'template_len_mean': 250, # mean length of template
  'template_len_sd': 30,    # sd of template length
  'max_p_error': 0.01,      # Maximum error rate at tip of read
  'k': 20,                  # exponent
  'gc_bias': {              # Omit if not modeling GC bias
    'bias_center': 0.5,     # Center of uni-modal window
    'bias_height': 1.5,     # Maximum departure from mean coverage (peak of curve)
    'bias_spread': 0.3      # Width measure of curve
  }
}
"""

_description = __doc__ + '\nExample parameter set:\n' + __example_param_text

_example_params = eval(__example_param_text)

from itertools import izip

import numpy as np

import mitty.lib.util as mutil
from mitty.plugins.reads.base_plugin import ReadModel

import logging
logger = logging.getLogger(__name__)


class Model(ReadModel):
  """Stock read plugin that approximates Illumina reads"""
  def __init__(self, read_len=100, template_len_mean=250, template_len_sd=50, max_p_error=0.01, k=20, gc_bias=None):
    """Initialization

    :param read_len:
    :param template_len_mean:
    :param template_len_sd:
    :param max_p_error:
    :param k:
    :param gc_bias: dictionary with keys
      bias_center' - Center of uni-modal window
      bias_height  - Maximum departure from mean coverage (peak of curve)
      bias_spread  - Width measure of curve
    :return:
    """
    self.read_len, self.template_len_mean, self.template_len_sd = read_len, template_len_mean, template_len_sd
    l = np.linspace(1e-10, 1, self.read_len)
    self.error_profile = max_p_error * (np.exp(k*l) - 1)/(np.exp(k) - 1)
    self.phred = ''.join([chr(int(33 + max(0, min(-10*np.log10(p), 93)))) for p in self.error_profile])
    self.gc_bias = gc_bias
    if gc_bias is not None:
      self.gc_curve = self.initialize_gc_curve()
    ReadModel.__init__(self, True)

  def initialize_gc_curve(self):
    """This creates a lookup table for the GC bias curve. We split the GC content range [0, 1.0] into 100 bins. We
    take the GC content of our read, multiply by 100, round to integer and look up the bias value (which is the
    modified probability of the read."""
    gc_f = np.linspace(0, 1.0, 100, dtype=float)
    center, height, spread = self.gc_bias['bias_center'], self.gc_bias['bias_height'], self.gc_bias['bias_spread']
    return height * np.exp(-((gc_f - center) / spread) ** 2)

  def gc_bias_reads(self, template_locs, template_lens, seq, gc_bias_rng):
    """

    :param template_locs:
    :param template_lens:
    :return:
    """
    rg = gc_bias_rng.rand(len(template_locs))
    s_cnt = seq.count
    gc_crv = self.gc_curve
    idx = [n for n, (r, t_loc, t_len) in enumerate(izip(rg, template_locs, template_lens))
           if r < gc_crv[int(100 * (s_cnt('G', t_loc, t_loc + t_len) + s_cnt('C', t_loc, t_loc + t_len)) / float(t_len))]]
    # return zip(*[(t_loc, t_len) for r, t_loc, t_len in izip(rg, template_locs, template_lens)
    #              if r < gc_crv[int(100 * (s_cnt('G', t_loc, t_loc + t_len) + s_cnt('C', t_loc, t_loc + t_len)) / float(t_len))]])
    return template_locs[idx], template_lens[idx]

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
    end_base = end_base or len(seq)
    #assert len(seq) > self.template_len_mean * 3, 'Template size should be less than 1/3rd sequence length'
    p_template = 0.5 * coverage / float(self.read_len)  # Per base probability of a template
    if self.gc_bias is not None:
      p_template *= self.gc_bias['bias_height']
    assert 0 < p_template < 0.1, 'Please lower coverage per block to less than {:f} otherwise accuracy will suffer'.format(0.1/p_template * coverage)

    template_loc_rng, read_order_rng, template_len_rng, error_loc_rng, base_choice_rng, gc_bias_rng = mutil.initialize_rngs(seed, 6)
    template_locs = mutil.place_poisson_seq(template_loc_rng, p_template, start_base, end_base, seq)
    template_lens = (template_len_rng.randn(template_locs.shape[0]) * self.template_len_sd + self.template_len_mean).astype('i4')
    idx = ((template_locs + template_lens < end_base) & (template_lens > self.read_len)).nonzero()[0]
    template_locs, template_lens = template_locs[idx], template_lens[idx]
    if self.gc_bias is not None:
      template_locs, template_lens = self.gc_bias_reads(template_locs, template_lens, seq, gc_bias_rng)

    read_order = read_order_rng.randint(2, size=template_locs.shape[0])  # Which read comes first?

    reads = np.recarray(dtype=ReadModel.dtype, shape=2 * template_locs.shape[0])
    r_start, r_o, pr = reads['start_a'], reads['read_order'], reads['perfect_reads']
    reads['read_len'] = self.read_len
    reads['read_order'][::2] = read_order[:]
    reads['read_order'][1::2] = 1 - read_order[:]
    r_len = self.read_len

    idx_fwd = (read_order == 0).nonzero()[0]
    idx_rev = read_order.nonzero()[0]

    r_start[2 * idx_fwd] = template_locs[idx_fwd]
    r_start[2 * idx_fwd + 1] = template_locs[idx_fwd] + template_lens[idx_fwd] - r_len
    r_start[2 * idx_rev] = template_locs[idx_rev] + template_lens[idx_rev] - r_len
    r_start[2 * idx_rev + 1] = template_locs[idx_rev]

    for n in xrange(reads.shape[0]):
      if r_o[n] == 0:  # Forward read comes first
        pr[n] = seq[r_start[n]:r_start[n] + r_len]
      else:  # Reverse strand read
        pr[n] = seq_c[r_start[n]:r_start[n] + r_len][::-1]

    if corrupt:
      self.corrupt_reads(reads, error_loc_rng, base_choice_rng)
    return reads, self.paired

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
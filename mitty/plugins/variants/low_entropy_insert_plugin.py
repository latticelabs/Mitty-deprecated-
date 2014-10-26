"""This generates insertions that consist of repeated copies of the same subsequence."""
import numpy
import mitty.plugins.variants.util as util
from mitty.lib.variation import Variation
import logging
logger = logging.getLogger(__name__)

__example_param_text = """
{
  "chromosome": [4],   # List of chromosomes to apply the variant to
  "p": 0.01,           # probability that the variant will occur at any given base
  "phet": 0.5,         # probability that the variant will be heterozygous
  "ins_len_lo": 10,
  "ins_len_hi": 20,
  "sub_seq_len": 5
}
"""

_description = """
This variant generates low complexity insertions that consist of repeated copies of the same subsequence.
""" + __example_param_text

_example_params = eval(__example_param_text)

#             0      1      2      3
gt_string = ['0/0', '0/1', '1/0', '1/1']  # The types of genotypes


def variant_generator(ref={},
                      chromosome=None,
                      p=0.01,
                      phet=0.5,
                      ins_len_lo=100,
                      ins_len_hi=10000,
                      sub_seq_len=10,
                      master_seed=1,
                      **kwargs):
  assert 0 <= p <= 1.0, "Probability out of range"
  assert 0 <= phet <= 1.0, "Probability out of range"
  logger.debug('Master seed: {:d}'.format(master_seed))
  ins_loc_rng, ins_len_rng, base_sel_rng, het_rng, copy_rng = util.initialize_rngs(master_seed, 5)

  for chrom in chromosome:
    ref_seq = ref[chrom]
    ins_locs, = numpy.nonzero(ins_loc_rng.rand(len(ref_seq)) < p)
    if ins_locs.size == 0: continue  # No variants here
    ins_lens = ins_len_rng.randint(low=ins_len_lo, high=ins_len_hi+1, size=ins_locs.size)
    het_type = util.het(ins_locs.size, phet, het_rng, copy_rng)

    yield {chrom:[Variation(pos + 1, pos + 2, ref_seq[pos], ref_seq[pos] +
           _repeat_sequence(ins_len, sub_seq_len, base_sel_rng), het)
           for het, pos, ins_len in zip(het_type, ins_locs, ins_lens) if ref_seq[pos] != 'N']}


def _repeat_sequence(seq_len=100, subseq_len=10, base_sel_rng=None, alphabet=['A', 'C', 'T', 'G']):
  """Create a sequence by repeating a sub-sequence

  .. note:: This is meant to be used internally

  Args:
    seq_len (int)         : Length of sequence.
    subseq_len (int)      : Length of repeating, random, sub-sequence
    base_sel_rng (object) : Random number generator e.g. numpy.random
    alphabet (list, str)  : (optional) List of characters constituting the alphabet

  Returns:
    str (str):  The sequence.
  """
  subseq = base_sel_rng.choice(alphabet, size=subseq_len, replace=True, p=[.3, .2, .2, .3]).tostring()
  subseq_len = len(subseq)
  return subseq * (seq_len / subseq_len) + subseq[:seq_len % subseq_len]


def test():
  """Basic test"""
  ref = {
    1: 'ACTGACTGACTG',
    2: 'ACTGACTGACTGACTGACTGACTG'
  }
  vg = variant_generator(ref, chromosome=[1, 2], p=1.0)
  for v_list in vg:
    assert isinstance(v_list.values()[0][0], Variation), v_list


if __name__ == "__main__":
  print _description
"""This generates insertions that consist of repeated copies of the same subsequence. The subsequence can be generated
de novo, or can be copied over from the neighborhood.

Example parameter snippet::

    "variant_models": [
      {
        "chromosome": [2],
        "model": "low_entropy_insert",
        "phet": 0.9,
        "p": 0.01,
        "ins_len_lo": 10,
        "ins_len_hi": 20,
        "sub_seq_len": 5,
        "copy_from_neighborhood": False,
        "master_seed": 0,
        "ins_loc_rng_seed": 1,
        "ins_len_rng_seed": 2,
        "base_sel_rng_seed": 3,
        "het_rng_seed": 4,
        "copy_rng_seed": 5
      },
      {
        "chromosome": [2],
        "model": "low_entropy_insert",
        "phet": 0.9,
        "p": 0.01,
        "ins_len_lo": 10,
        "ins_len_hi": 20,
        "sub_seq_len": 10,
        "copy_from_neighborhood": True,
        "master_seed": 0,
        "ins_loc_rng_seed": 1,
        "ins_len_rng_seed": 2,
        "base_sel_rng_seed": 3,
        "het_rng_seed": 4,
        "copy_rng_seed": 5
      }
    ]

"""
import numpy
import logging
import mitty.Plugins.Mutation.util as util
logger = logging.getLogger(__name__)
#             0      1      2      3
gt_string = ['0/0', '0/1', '1/0', '1/1']  # The types of genotypes


def _repeat_sequence(seq_len=100, subseq=None, subseq_len=10, base_sel_rng=None, alphabet=['A', 'C', 'T', 'G']):
  """Create a sequence by repeating a sub-sequence

  .. note:: This is meant to be used internally

  Parameters
  ----------
  seq_len      : int
                 Length of sequence.
  subseq       : str, optional
                 Sub-sequence to use as repeat block. Omit to generate a random sub-sequence.
  subseq_len   : int, optional
                 If subseq is omitted this must be provided to indicate desired length of random sub-sequence
  subseq_len   : int
                 Length of random sub-sequence
  base_sel_rng : object
                 Random number generator e.g. numpy.random
  alphabet     : list, optional
                 List of characters constituting the alphabet

  Returns
  -------
  str
      The sequence.

  Examples
  --------

  >>> _repeat_sequence(seq_len=7, subseq='ACTG')
  'ACTGACT'

  >>> _repeat_sequence(seq_len=7, subseq_len=4, base_sel_rng=numpy.random.RandomState(seed=1))
  'CGACCGA'

  """
  subseq = subseq or base_sel_rng.choice(alphabet, size=subseq_len, replace=True, p=[.3, .2, .2, .3]).tostring()
  subseq_len = len(subseq)
  return subseq * (seq_len / subseq_len) + subseq[:seq_len % subseq_len]


def _example_params():
  """This is used for the integration test."""
  return {
    "variant_models": [
      {
        "chromosome": [2],
        "model": "low_entropy_insert",
        "phet": 0.9,
        "p": 0.01,
        "ins_len_lo": 10,
        "ins_len_hi": 20,
        "sub_seq_len": 5,
        "copy_from_neighborhood": False,
        "master_seed": 0,
        "ins_loc_rng_seed": 1,
        "ins_len_rng_seed": 2,
        "base_sel_rng_seed": 3,
        "het_rng_seed": 4,
        "copy_rng_seed": 5
      },
      {
        "chromosome": [2],
        "model": "low_entropy_insert",
        "phet": 0.9,
        "p": 0.01,
        "ins_len_lo": 10,
        "ins_len_hi": 20,
        "sub_seq_len": 10,
        "copy_from_neighborhood": True,
        "master_seed": 0,
        "ins_loc_rng_seed": 1,
        "ins_len_rng_seed": 2,
        "base_sel_rng_seed": 3,
        "het_rng_seed": 4,
        "copy_rng_seed": 5
      }

    ]
  }


def variants(ref_fp=None,
             chromosome=None,
             p=0.01,
             phet=0.5,
             ins_len_lo=100,
             ins_len_hi=10000,
             sub_seq_len=10,
             copy_from_neighborhood=False,
             master_seed=None,
             ins_loc_rng_seed=1,
             ins_len_rng_seed=2,
             base_sel_rng_seed=3,
             het_rng_seed=4,
             copy_rng_seed=5,
             **kwargs):
  """Creates low entropy insertions.

  Parameters
  ----------
  ref_fp               : HDF5 file pointer
                         The reference whole genome file
  chromosome           : int list
                         All the chromosomes we want these variants to be sprinkled in
  p                    : float (0.0, 1.0)
                         Probability of inserts
  phet                 : float (0.0, 1.0)
                         Probability of having heterozygous mutation
  ins_len_lo           : int (> 0)
                         Shortest insertion length.
  ins_len_hi           : int (> ins_len_lo)
                         Longest insertion length.
  sub_seq_len          : int (> 0)
                         Length of the repeating block
  copy_from_neighborhood : bool
                         If True, instead of generating a random sequence that is repeated, a chunk of the reference
                         sequence neighboring the insert is used as the repeating block
  master_seed          : int
                         If this is specified then all individual seeds below are ignored and regenerated based on
                         `master_seed`. If this is not set (None) the individual seeds below are required
  ins_loc_rng_seed     : int
                         Seed for RNG determining insertion points
  ins_len_rng_seed     : int
                         Seed for RNG determining insertion lengths
  base_sel_rng_seed    : int
                         Seed for RNG determining bases in repeated sub-sequence.
                         Ignored if `copy_from_neighborhood` is True
  het_rng_seed         : int
                         Seed for RNG used to decide if genotype is heterozygous or not (0/1 or 1/0  vs 1/1)
  copy_rng_seed        : int
                         Seed for RNG used to decide which het type 0/1 or 1/0
  kwargs               : dict
                         Absorbs any other parameters it does not use

  Returns
  -------
  variant              : list of tuples
                         [(description, footprint, vcf_line) ... ]

  Notes
  -----
  - Insertion positions are picked via a poisson process.
  - Insertion lengths are uniformly distributed between lo and hi limits.
  """
  if master_seed:
    ins_loc_rng_seed, ins_len_rng_seed, base_sel_rng_seed, het_rng_seed, copy_rng_seed = \
      numpy.random.RandomState(seed=master_seed).randint(100000000, size=5)
    logger.debug('Used master seed to generate seeds {:d}, {:d}, {:d}, {:d}, {:d}'.
                 format(ins_loc_rng_seed, ins_len_rng_seed, base_sel_rng_seed, het_rng_seed, copy_rng_seed))
  ins_loc_rng, ins_len_rng, base_sel_rng, het_rng, copy_rng = \
    util.initialize_rngs(ins_loc_rng_seed, ins_len_rng_seed, base_sel_rng_seed, het_rng_seed, copy_rng_seed)

  description, footprint, vcf_line = [], [], []

  for chrom in chromosome:
    ref_seq = ref_fp['sequence/{:d}/1'.format(chrom)][:].tostring()  # Very cheap operation
    ins_locs, = numpy.nonzero(ins_loc_rng.rand(len(ref_seq)) < p)
    if ins_locs.size == 0: continue  # No variants here
    ins_lens = ins_len_rng.randint(low=ins_len_lo, high=ins_len_hi+1, size=ins_locs.size)
    het_type = util.het(ins_locs.size, phet, het_rng, copy_rng)

    for het, pos, ins_len in zip(het_type, ins_locs, ins_lens):
      ref = ref_seq[pos]
      if ref == 'N':
        continue  # Not valid, skip to next

      subseq = ref_seq[max(0, pos - sub_seq_len):pos] if copy_from_neighborhood else None
      # There is one edge case we are not handling here: if insert is located close to the beginning of the sequence
      # Too rare with real sequences
      alt = ref + _repeat_sequence(seq_len=ins_len, subseq=subseq, subseq_len=sub_seq_len, base_sel_rng=base_sel_rng)
      gt = gt_string[het]

      description.append('Insert')
      footprint.append([(het, chrom-1, pos, pos + 1)])  # Chrom is internal numbering, starts from 0
      # [(het, chrom, pos_st, pos_nd)]
      vcf_line.append([(chrom, pos+1, '.', ref, alt, 100, 'PASS', '.', 'GT', gt)])  # POS is VCF number starts from 1
      # CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample

    logger.debug('Generated {:d} insertions with min={:d},mean={:f},max={:d}'.
                 format(len(vcf_line), int(ins_lens.min()), ins_lens.mean(), int(ins_lens.max())))
  return description, footprint, vcf_line

if __name__ == "__main__":
  import pprint
  pprint.pprint(_example_params(), indent=1)
  print variants.__doc__
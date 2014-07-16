"""Given a whole genome, do calculations related to k-mers.

Usage:
  kmer  --wg=WG  --k=K  --prefout=OUT

Options:
  --wg=WG         Whole genome file name
  --k=K           k-mer length
  --prefout=OUT   Outfile name prefix
"""
import re
import docopt
import itertools
import numpy
import pylab
import cPickle
import genome
__version__ = '0.1.0'


def kmer_dict_for_seq(seq, kmer_len, alphabet='ACTG'):
  """Given a sequence and a kmer_length return a dictionary with all possible kmers and their counts.
  >>> kmer_dict_for_seq('ACTG', 1) == {'A': 1, 'C': 1, 'T': 1, 'G': 1}
  True
  >>> kd = kmer_dict_for_seq('ACTG', 2)
  >>> kd == {'AA': 0, 'CA': 0, 'TA': 0, 'GA': 0, 'AC': 1, 'CC': 0, 'TC': 0, 'GC': 0, 'AT': 0, 'CT': 1, 'TT': 0, 'GT': 0, 'AG': 0, 'CG': 0, 'TG': 1, 'GG': 0}
  True

  Any N's should be ignored
  >>> kd = kmer_dict_for_seq('ACTGNNNN', 2)
  >>> kd == {'AA': 0, 'CA': 0, 'TA': 0, 'GA': 0, 'AC': 1, 'CC': 0, 'TC': 0, 'GC': 0, 'AT': 0, 'CT': 1, 'TT': 0, 'GT': 0, 'AG': 0, 'CG': 0, 'TG': 1, 'GG': 0}
  True
  """
  return {''.join(k): seq.count(''.join(k)) for k in itertools.product(alphabet, repeat=kmer_len)}


def kmer_locate(seq, kmer_len, alphabet='ACTG'):
  """Given a sequence and a kmer_length return a dictionary with all possible kmers and their positions. We'd like to
  plot iki (inter-kmer-histograms) from this data.

  AGTCAGGNTCATGATTACTNTCAGTA
  01234567890123456789012345
  >>> kd = kmer_locate('AGTCAGGNTCATGATTACTNTCAGTA', 2, alphabet='ACTG')
  >>> kd == {'AA': [], 'CA': [3, 9, 21], 'TA': [15, 24], 'GA': [12], 'AC': [16], 'CC': [], 'TC': [2, 8, 20], 'GC': [], 'AT': [10, 13], 'CT': [17], 'TT': [14], 'GT': [1, 23], 'AG': [0, 4, 22], 'CG': [], 'TG': [11], 'GG': [5]}
  True

  AAANNAAAA
  012345678
  >>> kd = kmer_locate('AAANNAAAA', 2, alphabet='ACTG')
  >>> kd == {'AA': [0, 1, 5, 6, 7], 'CA': [], 'TA': [], 'GA': [], 'AC': [], 'CC': [], 'TC': [], 'GC': [], 'AT': [], 'CT': [], 'TT': [], 'GT': [], 'AG': [], 'CG': [], 'TG': [], 'GG': []}
  True

  """
  return {''.join(k): [m.start() for m in re.finditer('(?={:s})'.format(''.join(k)), seq)] for k in itertools.product(alphabet, repeat=kmer_len)}


def kmer_locate_sparse(seq, chrom_key, kmer_len, kmer_dict=None, alphabet='ACTG'):
  """Given a sequence and a kmer_length return a dictionary with only existing kmers and their positions.

  AAATTTCCCGGGAAATTTCCCGGG
  012345678901234567890123
  >>> kd = kmer_locate_sparse('AAATTTCCCGGGAAATTTCCCGGG', chrom_key='1:1', kmer_len=2, alphabet='ACTG')
  >>> kd == {'AA': {'1:1': [0, 1, 12, 13]}, 'AT': {'1:1': [2, 14]},
  ... 'TT': {'1:1': [3, 4, 15, 16]}, 'TC': {'1:1': [5, 17]}, 'CC': {'1:1': [6, 7, 18, 19]},
  ... 'CG': {'1:1': [8, 20]}, 'GG': {'1:1': [9, 10, 21, 22]}, 'GA': {'1:1': [11]}}
  True

  ATCGATCGATCG
  012345678901
  >>> kd = kmer_locate_sparse('ATCGATCGATCG', chrom_key='2:1', kmer_len=2, alphabet='ACTG', kmer_dict=kd)
  >>> kd == {'AA': {'1:1': [0, 1, 12, 13]}, 'CC': {'1:1': [6, 7, 18, 19]}, 'TT': {'1:1': [3, 4, 15, 16]},
  ... 'CG': {'1:1': [8, 20], '2:1': [2, 6, 10]}, 'GG': {'1:1': [9, 10, 21, 22]}, 'AT': {'1:1': [2, 14], '2:1': [0, 4, 8]},
  ... 'GA': {'1:1': [11], '2:1': [3, 7]}, 'TC': {'1:1': [5, 17], '2:1': [1, 5, 9]}}
  True
  """
  if kmer_dict is None:
    kmer_dict = {}
  for n in range(len(seq) - kmer_len + 1):
    key = seq[n:n + kmer_len]
    if 'N' in key:
      continue
    try:
      k = kmer_dict[key]
    except KeyError:
      k = kmer_dict[key] = {}
    try:
      k[chrom_key].append(n)
    except KeyError:
      k[chrom_key] = [n]

  return kmer_dict


# def kmer_candidate_loc_count(wg, chrom_key, kmer_len):
#   """Given a sequence dictionary (like a WholeGenome) return us a list of locations along with counts of how many
#   locations in the genome the k-mer starting there maps to as well. Minimum should be 1. If it is 0 it means that the
#   sequence there contains an 'N' and was discarded
#   """
#   this_seq = wg[chrom_key][0]
#   candidate_count = numpy.recarray((len(this_seq) - kmer_len + 1,), dtype=[(k,int) for k in wg.index.keys()])
#   for n in range(len(this_seq) - kmer_len + 1):  # Walk along this sequence
#     this_kmer = this_seq[n:n + kmer_len]
#     for chrom in wg.index.keys():
#       candidate_count[chrom][n] = 0 if 'N' in this_kmer else wg[chrom][0].count(this_kmer)
#
#   return candidate_count

def kmer_candidate_loc_count(wg, chrom_key, kmer_len):
  """Given a sequence dictionary (like a WholeGenome) return us a list of locations along with counts of how many
  locations in the genome the k-mer starting there maps to as well. Minimum should be 1. If it is 0 it means that the
  sequence there contains an 'N' and was discarded
  """
  this_seq = wg[chrom_key][0]
  candidate_count = numpy.zeros(len(this_seq) - kmer_len + 1, dtype=numpy.uint32)  # We are being safe
  for n in range(len(this_seq) - kmer_len + 1):  # Walk along this sequence
    this_kmer = this_seq[n:n + kmer_len]
    if 'N' in this_kmer:
      continue
    for chrom in wg.index.keys():
      candidate_count[n] += wg[chrom][0].count(this_kmer)
  return candidate_count


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  #Load the whole genome file
  with genome.WholeGenome(fname=args['--wg']) as wg:
    for chrom_key in wg.index.keys():
      d = kmer_candidate_loc_count(wg, chrom_key=chrom_key, kmer_len=int(args['--k']))
      with open(args['--prefout'] + '_' + chrom_key + '.pkl', 'w') as fout:
        cPickle.dump(d, fout, protocol=cPickle.HIGHEST_PROTOCOL)
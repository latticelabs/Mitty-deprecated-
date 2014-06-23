"""Given a whole genome, do calculations related to k-mers.

Usage:
  kmer  --wg=WG  --k=K  --out=OUT

Options:
  --wg=WG    Whole genome file name
  --k=K       k-mer length
  --out=OUT  Outfile name
"""
import re
import docopt
import itertools
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


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  #Load the whole genome file
  with genome.WholeGenome(fname=args['--wg']) as wg:
    k_mer_locations = {k: kmer_locate(seq=wg[k][0], kmer_len=int(args['--k'])) for k in wg.index.keys()}
    with open(args['--out'], 'w') as fout:
      cPickle.dump(k_mer_locations, fout, protocol=cPickle.HIGHEST_PROTOCOL)
"""Given a whole genome, do calculations related to k-mers.

Usage:
  kmer  --wg=WG  --k=K  --out=OUT

Options:
  --wg=WG    Whole genome file name
  --k=K       k-mer length
  --out=OUT  Outfile name
"""
import docopt
import cPickle
import genome
__version__ = '0.1.0'


def kmer_iter(length, alphabet=['A', 'C', 'T', 'G']):
  """
  >>> my_kmer = kmer_iter(length=2)
  >>> [k for k in my_kmer]
  ['AA', 'CA', 'TA', 'GA', 'AC', 'CC', 'TC', 'GC', 'AT', 'CT', 'TT', 'GT', 'AG', 'CG', 'TG', 'GG']
  """
  l_alph = len(alphabet)
  kmer_count = l_alph ** length
  cntr = [-1, 0, 0, 0]
  for n in range(kmer_count):
    digit_index = 0
    while True:
      cntr[digit_index] += 1
      if cntr[digit_index] == l_alph:
        cntr[digit_index] = 0
        digit_index += 1
      else:
        break
    yield ''.join([alphabet[cntr[dig]] for dig in range(length)])


def kmer_dict_for_seq(seq, kmer_len):
  """Given a sequence and a kmer_length return a dictionary with all possible kmers and their counts.
  >>> kmer_dict_for_seq('ACTG', 1)
  {'A': 1, 'C': 1, 'T': 1, 'G': 1}
  >>> kd = kmer_dict_for_seq('ACTG', 2)
  >>> kd == {'AA': 0, 'CA': 0, 'TA': 0, 'GA': 0, 'AC': 1, 'CC': 0, 'TC': 0, 'GC': 0, 'AT': 0, 'CT': 1, 'TT': 0, 'GT': 0, 'AG': 0, 'CG': 0, 'TG': 1, 'GG': 0}
  True

  Any N's should be ignored
  >>> kd = kmer_dict_for_seq('ACTGNNNN', 2)
  >>> kd == {'AA': 0, 'CA': 0, 'TA': 0, 'GA': 0, 'AC': 1, 'CC': 0, 'TC': 0, 'GC': 0, 'AT': 0, 'CT': 1, 'TT': 0, 'GT': 0, 'AG': 0, 'CG': 0, 'TG': 1, 'GG': 0}
  True
  """
  my_kmer = kmer_iter(length=kmer_len)
  kmer_dict = {k: 0 for k in my_kmer}
  for n in range(len(seq) - kmer_len + 1):
    try:
      kmer_dict[seq[n:n + kmer_len]] += 1
    except KeyError:  # This will be because we have N's some where. Ignore this stuff
      continue
  return kmer_dict

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  #Load the whole genome file
  with genome.WholeGenome(fname=args['--wg']) as wg:
    k_mer_counts = {k: None for k in wg.index.keys()}
    for chrom in wg.index.keys():
      print chrom
      k_mer_counts[chrom] = kmer_dict_for_seq(seq=wg[chrom][0], kmer_len=int(args['--k']))
    with open(args['--out'], 'w') as fout:
      cPickle.dump(k_mer_counts, fout, protocol=cPickle.HIGHEST_PROTOCOL)
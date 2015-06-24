"""
Command-line program and library for counting k-mers in a genome.

Produces a Python pickle file with a dictionary of k-mer counts in the genome::
  kmers = {
    'md5': {
      'gi|568336015|gb|CM000671.2| Homo sapiens chromosome 9, GRCh38 reference primary assembly': '6c198acf68b5af7b9d676dfdd531b5de',
      'gi|568336001|gb|CM000685.2| Homo sapiens chromosome X, GRCh38 reference primary assembly': '2b3a55ff7f58eb308420c8a9b11cac50',
      ...
    }
    'kmers': {
      'GATGACATGG': 2993,
      'TTATTTATTG': 19149,
      ...
    }
  }
(The 'md5' entry stores the sequence ids and md5 sums so we can cross check with the reference genome if we wish.)

If asked to compute a genome score, a sampled k-mer score is computed for each chromosome and is dumped as a list of
arrays (one for every chromosome) in a .npz (compressed numpy) file. The first 'array' in the list has a single element
corresponding to the value of 'D' (score computed every D bases)

Commandline::

  Usage:
    kmers count <fasta> <kfile> [-k=K] [-v] [-p]
    kmers score <fasta> <kfile> <sfile> [-d=D] [-v] [-p]

  Options:
    count     Create a table of k-mer counts in the genome
    score     Score whole genome using k-mer counts
    <fasta>   Fasta file
    <kfile>   Name of the k-mer count table file (Python .pkl file).
    -k=K      k-mer size [default: 10]
    <sfile>   Genome score file (If omitted this will not be computed)
    -d=D      Decimate bases. Compute k-score every D bases [default: 1000]
"""
import time
import cPickle

import docopt
import numpy as np

import mitty.lib.io as mio
import mitty.lib
import mitty.lib.util as mutil

import logging
logger = logging.getLogger(__name__)


def process_genome(ref, k=10, progress_bar_func=None):
  """Compute kmer dictionary given a genome and value for k. Also give us the md5sums."""
  kmers = {}
  md5 = {}
  if progress_bar_func: progress_bar_func('Counting ', 0, 40)
  for n, seq in enumerate(ref):
    md5[seq['id']] = seq['md5']
    mutil.parse_sequence(seq['seq'], k=10, kmers=kmers)
    if progress_bar_func: progress_bar_func('Counting ', (n + 1) / float(len(ref)), 40)
  if progress_bar_func:
    progress_bar_func('Counting ', 1.0, 40)
    print('')
  return kmers, md5


def score_entire_genome(ref, k_mer_table, d=100, progress_bar_func=None):
  """Given a sequence compute the k-mer scores at selected points

  :param ref:
  :param k_mer_table:
  :param d:  Bases to skip when computing k-score
  :return: A list of numpy arrays
  """
  k = len(k_mer_table.keys()[0])
  scores = {}
  for n, seq in enumerate(ref):
    if progress_bar_func: progress_bar_func('Scoring ', (n + 1) / float(len(ref)), 40)
    scores['chrom_{:d}'.format(n + 1)] = mutil.score_long_sequence(seq['seq'], d, k, k_mer_table)
  if progress_bar_func:
    progress_bar_func('Scoring ', 1.0, 40)
    print('')
  return scores


def cli():
  """Serves as entry point for scripts"""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__)  # Version string?

  logging.basicConfig()
  if cmd_args['-v']:
    logger.setLevel(logging.DEBUG)
  progress_bar_func=mitty.lib.progress_bar if cmd_args['-p'] else None
  ref = mio.Fasta(multi_fasta=cmd_args['<fasta>'])
  k_f_name = cmd_args['<kfile>']
  k = int(cmd_args['-k'])
  s_f_name = cmd_args['<sfile>']
  d = int(cmd_args['-d'])

  if cmd_args['count']:  # Compute k-mer table
    count(ref, k_f_name=k_f_name, k=k, progress_bar_func=progress_bar_func)
  elif cmd_args['score']:  # Score the genome
    score(ref, k_f_name=k_f_name, s_f_name=s_f_name, d=d, progress_bar_func=progress_bar_func)


def count(ref, k_f_name, k, progress_bar_func=None):
  t0 = time.time()
  k_mer_count_table, md5 = process_genome(ref, k=k, progress_bar_func=progress_bar_func)
  t1 = time.time()
  logger.debug('Counting: {:f} s'.format(t1 - t0))
  logger.debug('{:d}-mer coverage: {:d}/{:d}'.format(k, len(k_mer_count_table.keys()), 4 ** k))
  with open(k_f_name, 'wb') as fp:
    cPickle.dump({'md5': md5, 'kmers': k_mer_count_table}, fp, protocol=cPickle.HIGHEST_PROTOCOL)


def score(ref, k_f_name, s_f_name, d, progress_bar_func=None):
  saved_k_mer_table_data = cPickle.load(open(k_f_name, 'rb'))
  k_mer_count_table, md5 = saved_k_mer_table_data['kmers'], saved_k_mer_table_data['md5']
  k = len(k_mer_count_table.keys()[0])
  logger.debug('{:d}-mer coverage: {:d}/{:d}'.format(k, len(k_mer_count_table.keys()), 4 ** k))
  t0 = time.time()
  scores = score_entire_genome(ref, k_mer_count_table, d=d, progress_bar_func=progress_bar_func)
  scores['step'] = [d]
  np.savez_compressed(s_f_name, **scores)
  t1 = time.time()
  logger.debug('Scoring: {:f}'.format(t1 - t0))


if __name__ == '__main__':
  cli()
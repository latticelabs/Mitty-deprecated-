"""Command-line program and library for counting k-mers in a genome. Produces a Python pickle file with a dictionary
of k-mer counts in the genome. The 'md5' entry stores the sequence ids and md5 sums so we can cross check with the
reference genome if we wish.
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

Commandline::

  Usage:
    countkmers  <fasta> <out> [-k=K] [-v] [-p]

  Options:
    <fasta>     Fasta file
    -k=K        k-mer size [default: 10]
    <out>       Output file name (Will be a Python pkl file)
    -p          Show progress bar
    -v          Dump detailed logger messages

"""
import time
import cPickle

import docopt

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


def average_kmer_counts(seq_list, kmers):
  """Given a list of sequences given a list of values indicating the average k-mer count of that sequence"""
  if len(kmers) == 0: return None
  avg_counts = []
  k = len(kmers.keys()[0])
  for seq in seq_list:
    score = 0
    segments = 0
    for n in xrange(0, len(seq) - k):
      score += kmers[seq[n:n+k]]
      segments += 1
    avg_counts.append(score/float(segments))
  return avg_counts


def cli():
  """Serves as entry point for scripts"""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__)  # Version string?

  logging.basicConfig()
  if cmd_args['-v']:
    logger.setLevel(logging.DEBUG)

  ref = mio.Fasta(multi_fasta=cmd_args['<fasta>'])
  k = int(cmd_args['-k'])
  t0 = time.time()
  kmers, md5 = process_genome(ref, k=k, progress_bar_func=mitty.lib.progress_bar if cmd_args['-p'] else None)
  t1 = time.time()
  logger.debug('{:f} s'.format(t1 - t0))

  with open(cmd_args['<out>'], 'wb') as fp:
    cPickle.dump({'md5': md5, 'kmers': kmers}, fp, protocol=cPickle.HIGHEST_PROTOCOL)

  logger.debug('{:d}-mer coverage: {:d}/{:d}'.format(k, len(kmers.keys()), 4 ** k))

if __name__ == '__main__':
  cli()
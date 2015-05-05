"""Command-line program and library for counting k-mers in a genome.

Commandline::

  Usage:
    kmer  <fasta> <out> [-k=K] [-v] [-p]

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

  t0 = time.time()
  kmers, md5 = process_genome(ref, k=int(cmd_args['-k']),
                              progress_bar_func=mitty.lib.progress_bar if cmd_args['-p'] else None)
  t1 = time.time()
  logger.debug('{:f} s'.format(t1 - t0))

  with open(cmd_args['<out>'], 'wb') as fp:
    cPickle.dump({'md5': md5, 'kmers': kmers}, fp, protocol=cPickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
  cli()
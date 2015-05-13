#!python
"""Given a gzipped multi fa.gz file split it into separate files as used by Mitty

Commandline::

  Usage:
    splitta  <fagz>  <dout>  [-v]

  Options:
    <fagz>   Fasta Gzip file input
    <dout>   Output directory
    -v       Dump detailed logger messages
"""
__version__ = '1.0.0'

import os
import hashlib

import docopt

import logging
logger = logging.getLogger(__name__)


def split_multi_fasta_gz(fa_fname, dir_out):
  """Given a gzipped multi fa.gz file split it into separate files as used by mitty.
  Index file contains file seqid, len and md5 sum in same order as in fasta.gz file"""
  def write_it_out(d_out, ch, sid, seq):
    logger.debug('Writing out {:s}'.format(sid))
    with open(os.path.join(d_out, 'chr{:d}.fa'.format(ch)), 'w') as fp_out:
      fp_out.write(sid + '\n')
      s = ''.join(seq)
      fp_out.write(s)
    with open(os.path.join(d_out, 'index.csv'), 'a') as fp_out:
      fp_out.write('{:s}\t{:d}\t{:s}\n'.format(sid, len(s), hashlib.md5(s).hexdigest()))

  import gzip
  import os
  with gzip.open(fa_fname, 'rb') as fp:
    this_seq = []
    chrom = 0
    try:
      os.remove(os.path.join(dir_out, 'index.csv'))
    except OSError:
      pass
    for line in fp:
      line = line.strip()
      if len(line) == 0: continue
      if line[0] == '>':
        if len(this_seq):
          write_it_out(dir_out, chrom, seq_id, this_seq)
        seq_id = line[1:]
        this_seq = []
        chrom += 1
      else:
        this_seq += [line.upper()]
    if len(this_seq):
      write_it_out(dir_out, chrom, seq_id, this_seq)


def cli():
  """Serves as entry point for scripts"""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  if not os.path.exists(cmd_args['<dout>']):
    os.makedirs(cmd_args['<dout>'])

  split_multi_fasta_gz(cmd_args['<fagz>'], cmd_args['<dout>'])

if __name__ == "__main__":
  cli()
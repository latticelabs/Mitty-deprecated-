#!python
"""Given a gzipped multi fa.gz file split it into separate files as used by Mitty

Commandline::

  Usage:
    splitta  --fagz=FAGZ  --dir_out=DIR_OUT  [-v]

  Options:
    --fagz=FAGZ             Fasta Gzip file input
    --dir_out=DIR_OUT       Output directory
    -v                      Dump detailed logger messages
"""
__version__ = '1.0.0'

import os
import docopt
import logging
logger = logging.getLogger(__name__)


def split_multi_fasta_gz(fa_fname, dir_out):
  """Given a gzipped multi fa.gz file split it into separate files as used by mitty"""
  def write_it_out(d_out, ch, sid, seq):
    logger.debug('Writing out {:s}'.format(sid))
    with open(os.path.join(d_out, 'chr{:d}.fa'.format(ch)), 'w') as fp_out:
      fp_out.write(sid + '\n')
      fp_out.writelines(seq)

  import gzip
  import os
  with gzip.open(fa_fname, 'rb') as fp:
    this_seq = []
    chrom = 0
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

  if not os.path.exists(cmd_args['--dir_out']):
    os.makedirs(cmd_args['--dir_out'])

  split_multi_fasta_gz(cmd_args['--fagz'], cmd_args['--dir_out'])

if __name__ == "__main__":
  cli()
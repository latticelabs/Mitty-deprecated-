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
from mitty.lib.genome import split_multi_fasta_gz
import logging
logger = logging.getLogger(__name__)


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
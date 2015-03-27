#!python
"""Inspect genome databases created by `genomes` program

Commandline::

  Usage:
    inspect --dbfile=DBFILE
    inspect sfs <chrom>  --dbfile=DBFILE

  Options:
    sfs     Print the actual site frequency spectrum of the master list
    chrom   Which chrom to inspect
"""
import docopt
import numpy as np

import mitty.lib.db as mdb
import mitty.lib.variants as vr


def cli():  # pragma: no cover
  """Serves as entry point for scripts"""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__)  # Version string?

  if cmd_args['sfs']:
    print_sfs(cmd_args)
  else:
    inspect(cmd_args)


def inspect(cmd_args):
  """Print some useful information about the database

  :param cmd_args: parsed arguments
  """
  pop_db_name = cmd_args['--dbfile']
  conn = mdb.connect(db_name=pop_db_name)
  chrom = mdb.chromosomes_in_db(conn)
  n_s = mdb.samples_in_db(conn)
  variant_stats = np.empty((len(chrom), 3), dtype=float)
  sample_max = 100
  for i, c in enumerate(chrom):
    if n_s < sample_max:  # Take every sample
      ss = range(n_s)
    else:
      ss = np.random.randint(0, n_s, sample_max)
    s_len = np.empty(len(ss), dtype=float)
    for j, s in enumerate(ss):
      s_len[j] = len(mdb.load_sample(conn, 0, s, c))
    _, seq_len = mdb.load_chromosome_metadata(conn, c)
    variant_stats[i, :] = (s_len.mean(), s_len.std(), 1e6 * s_len.mean() / float(seq_len))

  print('{:s}'.format(pop_db_name))
  print('\t{:d} chromosomes'.format(len(chrom)))
  print('\t{:d} samples'.format(n_s))
  print('Population')
  print('\tChrom\tVariants')
  for c in chrom:
    print('\t{:d}\t{:d}'.format(c, mdb.variants_in_master_list(conn, c)))
  print('Samples')
  print('\tChrom\tAvg variants\tStd variants\tVariants/megabase')
  for i, c in enumerate(chrom):
    print('\t{:d}\t{:<9.2f}\t{:<9.2f}\t{:.1f}'.format(c, variant_stats[i, 0], variant_stats[i, 1], variant_stats[i, 2]))


def print_sfs(cmd_args):
  pop_db_name = cmd_args['--dbfile']
  conn = mdb.connect(db_name=pop_db_name)
  chrom = int(cmd_args['<chrom>'])
  ml = mdb.load_master_list(conn, chrom)
  print('Site frequency spectrum for chrom {:d}'.format(chrom))
  print(ml)
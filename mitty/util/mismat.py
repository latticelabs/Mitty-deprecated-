"""Program that takes one or more mis-alignment databases (as produced by perfectbam) and produces a summary matrix
comparing the mis-alignments in the files.

Commandline::

  Usage:
    mismat <misdb>...  [-c=C] [--vars] [-r=R] [--html]

  Options:
    <out>       Name of output file
    <misdb>     List of mis-alignment database files we want to summarize
    -c=C        Chrom number. Leave out to do across all chroms
    --vars      Set this option to compute data from reads containing variants only
    -r=R        Get reads from cell r,c,i r = [0,1,2 ..] c = [0, 1,2 ..] i = 'i' or 'd'
                There can be no spaces in the argument, or enclose in quotes
                e.q. -r 1,2,i  OR -r "1, 2, i"
    --html      Set this option to dump output as interactive html file
"""
import docopt
import sqlite3 as sq

import mitty.util.mis_db as mdb


def cli():
  """Command line entry point"""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
    return
  else:
    args = docopt.docopt(__doc__)

  conn = sq.connect(':memory:')
  mdb.attach_databases(conn, args['<misdb>'])
  db_count = len(args['<misdb>'])
  chrom = int(args['-c']) if args['-c'] else None
  variants_only = args['--vars']
  bm_mat = mdb.get_bench_mat(conn, db_count, chrom=chrom, variants_only=variants_only)
  print_benchmark_mat(bm_mat)
  if args['-r']:
    r, c, i = args['-r'].split(',')
    r, c, i = int(r), int(c), i.strip()
    assert 0 <= r < db_count
    assert 0 <= c < db_count
    assert i in ['i', 'd']
    for read in mdb.get_reads_mat_cell(conn, r, c, i, chrom=chrom, variants_only=variants_only):
      print(read)


def print_benchmark_mat(bm_mat):
  db_count = len(bm_mat['metadata']['dbs'])
  print('')
  for c in range(db_count):
    print(str(c) + ' = ' + bm_mat['metadata']['dbs']['db' + str(c)])
  if bm_mat['metadata']['chrom']: print('chrom = ' + str(bm_mat['metadata']['chrom']))
  print('Key: diff (row - col) / intersection (row ^ col)')

  print('\t' + ''.join(['\t' * 2 * c + '|' + str(c) for c in range(db_count)]))
  print('-' * 8 + '-' * 16 * db_count)
  for r in range(db_count):
    print(str(r) + '\t' + ''.join(['\t' * 2 * c + '|' + str(bm_mat.get((r, c, 'd'), '-')) + '/' + str(bm_mat.get((r, c, 'i'), '-')) for c in range(db_count)]))


if __name__ == '__main__':
  cli()
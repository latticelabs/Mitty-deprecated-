#!python
__cmd__ = """Commandline::

  Usage:
    perfectbam <inbam> <outbam> <dbname> [--window=WN] [-x] [-v] [-p]

  Options:
    <inbam>        Input bam file name of reads
    <outbam>       Perfect BAM will be written to this file
    <dbname>       Name of database file to store results in
    --window=WN    Size of tolerance window [default: 0]
    -x             Use extended CIGAR ('X's and '='s) rather than traditional CIGAR (just 'M's)
    -v             Dump detailed logger messages
    -p             Show progress bar
"""
__param__ = """Given a bam file containing simulated reads aligned by a tool:
  1. Produce a new bam that re-aligns all reads so that their alignment is perfect
  2. Produce a database file containing data about the misaligned and unmapped reads with the following columns::
       qname              -  qname of the read
       error_type         -  type of error 3 bit number  bit 0=chrom, 1=pos, 2=cigar
       correct_chrom      -  correct chromosome number of read
       correct_pos        -  correct position of read
       correct_cigar      -  correct cigar of read
       aligned_chrom      -  actual aligned chromosome
       aligned_pos        -  actual aligned pos
       aligned_cigar      -  actual aligned cigar
       mapping_qual       -  mapping quality
       mate_is_unmapped   -  is mate unmapped
       seq                -  actual sequence string"""
__doc__ = __cmd__ + __param__

import os
import sqlite3 as sq
import time

import pysam
from collections import Counter
import docopt

import mitty.lib.io as mio  # For the bam sort and index function
from mitty.lib.reads import old_style_cigar
from mitty.lib import progress_bar

import logging
logger = logging.getLogger(__name__)


def connect_to_db(db_name):
  """Convenience function so we don't need to import sqlite in other modules just for this"""
  conn = sq.connect(db_name)
  conn.row_factory = sq.Row
  return conn


def create_db(conn):
  """Create tables for the mis-aligned reads database

  :param conn: database connection
  """
  c = conn.cursor()
  c.execute('CREATE TABLE reads (qname TEXT, error_type INT, '
            'correct_chrom INT, correct_pos INT, correct_cigar TEXT, '
            'aligned_chrom INT, aligned_pos INT, aligned_cigar TEXT,'
            'mapping_qual INT, mate_is_unmapped BOOL, seq TEXT)')
  c.execute('CREATE TABLE summary (chrom INT, total_reads INT, incorrect_reads INT, unmapped_reads INT, seq_len INT, seq_id TEXT)')
  conn.commit()


def write_summary_to_db(conn, total_reads_cntr, incorrectly_aligned_reads_cntr, unmapped_reads_cntr, bam_seq_header):
  """Write the alignment analysis summary for the BAM file into the database

  :param conn: database connection
  :param total_reads_cntr: Counter by chromosome
  :param incorrectly_aligned_reads_cntr: Counter by chromosome
  :param unmapped_reads_cntr: Counter by chromosome
  :param bam_seq_header: bam_in_fp.header['SQ']
  :return:
  """
  insert_clause = "INSERT INTO summary VALUES (?, ?, ?, ?, ?, ?)"
  to_insert = [(n + 1, total_reads_cntr[n + 1], incorrectly_aligned_reads_cntr[n + 1], unmapped_reads_cntr[n + 1], seq['LN'], seq['SN']) for n, seq in enumerate(bam_seq_header)]
  conn.executemany(insert_clause, to_insert)


def write_reads_to_db(conn, data_to_save):
  """Write mis-aligned/unmapped read data to the database

  :param conn: db connection
  :param data_to_save: [qname, error_type, chrom, pos, cigar, aligned_chrom, aligned_pos, aligned_cigar,
                        map_quality, mate_is_unmapped, sequence]
  :return:
  """
  insert_clause = "INSERT INTO reads VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
  conn.execute(insert_clause, data_to_save)


def commit_and_create_db_indexes(conn):
  """Need to commit the INSERTs and create search indexes for faster querying

  :param conn: db connection
  """
  conn.commit()
  conn.execute('CREATE INDEX idx_f ON reads (correct_chrom, correct_pos)')
  conn.execute('CREATE INDEX idx_r ON reads (aligned_chrom, aligned_pos)')


def main(bam_in_fp, bam_out_fp, db_name, window, extended=False, progress_bar_func=None):
  """Main processing function that goes through the bam file, analyzing read alignment and writing out

  :param bam_in_fp:  Pointer to original BAM
  :param bam_out_fp: Pointer to perfect BAM being created
  :param db_name:    database name
  :param window:     Tolerance window for deciding if read is correctly aligned
  :param extended:   If True write out new style CIGARs (With '=' and 'X')
  :param progress_bar_func: Our famous progress bar, if we want it
  :return: number of reads processed
  """
  conn = sq.connect(db_name)
  create_db(conn)

  total_read_count = float(bam_in_fp.mapped + bam_in_fp.unmapped)
  f0 = 0
  total_reads_cntr, incorrectly_aligned_reads_cntr, unmapped_reads_cntr = Counter(), Counter(), Counter()

  if progress_bar_func is not None: progress_bar_func('Processing BAM ', 0, 80)

  for n, read in enumerate(bam_in_fp):
    # qname = '{:d}|{:d}|{:d}|{:d}|{:s}|{:d}|{:s}'
    try:
      if read.is_paired:
        if read.is_read1:
          _, chrom, cpy, ro, pos, cigar, _, _, _ = read.qname.split('|')
        else:
          _, chrom, cpy, _, _, _, ro, pos, cigar = read.qname.split('|')
      else:
        _, chrom, cpy, ro, pos, cigar = read.qname.split('|')[:6]  # For Wan-Ping :)
      ro, chrom, pos = int(ro), int(chrom), int(pos)
    except ValueError:
      logger.debug('Error processing qname: n={:d}, qname={:s}, chrom={:d}, pos={:d}'.format(n, read.qname, read.reference_id + 1, read.pos))
      continue

    if not extended:
      cigar = old_style_cigar(cigar)

    total_reads_cntr[chrom] += 1

    error_type = 0x0
    if read.is_unmapped:
      error_type = 0x8
    else:
      if read.reference_id != chrom - 1:
        error_type |= 0x1
      if not (-window <= read.pos - pos <= window):
        error_type |= 0x2
      if read.cigarstring != cigar:
        error_type |= 0x4
    if error_type != 0:
      if error_type == 0x8:  # Unmapped read
        unmapped_reads_cntr[chrom] += 1
      else:
        incorrectly_aligned_reads_cntr[chrom] += 1
      write_reads_to_db(conn,
                        [read.qname, error_type, chrom, pos, cigar,
                         read.reference_id + 1, read.pos, read.cigarstring,
                         read.mapq, read.mate_is_unmapped, read.query_sequence])

    # Now write out the perfect alignment
    read.is_reverse = 1 - ro
    read.mate_is_reverse = ro
    read.mate_is_unmapped = False  # Gotta check this - what if mate is deep in an insert?
    read.reference_id = chrom - 1
    read.pos = pos
    read.cigarstring = cigar  # What if this is deep in an insert?
    bam_out_fp.write(read)

    if progress_bar_func is not None:
      f = n / total_read_count
      if f - f0 >= 0.01:
        progress_bar_func('Processing BAM ', f, 80)
        f0 = f
  if progress_bar_func is not None: print('\n')

  write_summary_to_db(conn, total_reads_cntr, incorrectly_aligned_reads_cntr, unmapped_reads_cntr, bam_in_fp.header['SQ'])
  commit_and_create_db_indexes(conn)
  return int(total_read_count)


def load_summary(conn):
  """Load information in the summary table

  :param conn: database connection object
  :return:
  """
  return [{k: row[k] for k in row.keys()} for row in conn.execute('SELECT * FROM summary ORDER BY chrom')]


def load_reads(conn, chrom=None, start_pos=None, stop_pos=None, source=True, sample=None):
  """Load mis-aligned reads from the database based on the filters we give it

  :param conn: db connection object
  :param chrom: chromosome number [1,2...]. If None, load all chroms and ignore start_pos and stop_pos values
  :param start_pos: If chrom is given, start taking reads from this here
  :param stop_pos: If chrom given, stop reads here
  :param source: If True find reads with correct position matching criteria. If False find reads with aligned pos matching criteria
  :param sample: If given, sub-sample reads to get at most this many
  :return: list of tuples (corr_chrom, corr_pos, align_chrom, align_pos)

  SELECT * FROM reads WHERE rowid in (SELECT abs(random()) % (select max(rowid) FROM reads) FROM reads LIMIT 10);
  """
  select_statement = 'SELECT * FROM reads'
  col_prefix = 'correct_' if source else 'aligned_'
  where_clause = []
  if chrom is not None:
    where_clause += [col_prefix + 'chrom = {:d}'.format(chrom)]
    if start_pos is not None:
      where_clause += [col_prefix + 'pos >= {:d} '.format(start_pos)]
    if stop_pos is not None:
      where_clause += [col_prefix + 'pos <= {:d} '.format(stop_pos)]
  if sample is not None:
    where_clause += ['rowid in (SELECT abs(random()) % (select MAX(rowid) FROM reads) FROM reads LIMIT {:d})'.format(sample)]
  query = select_statement + (' WHERE ' if where_clause else '') + ' AND '.join(where_clause)
  return [r for r in conn.execute(query)]


def cli():
  """Command line script entry point."""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__cmd__, ['-h'])
  else:
    args = docopt.docopt(__cmd__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  with pysam.AlignmentFile(args['<inbam>'], 'rb') as bam_in_fp, \
      pysam.AlignmentFile(args['<outbam>'], 'wb', template=bam_in_fp) as bam_out_fp:
    try:
      os.remove(args['<dbname>'])
    except OSError:
      pass
    t0 = time.time()
    read_count = main(bam_in_fp=bam_in_fp, bam_out_fp=bam_out_fp, db_name=args['<dbname>'],
                      window=int(args['--window']), extended=bool(args['-x']), progress_bar_func=progress_bar if args['-p'] else None)
    t1 = time.time()
    logger.debug('Analyzed {:d} reads in BAM in {:2.2f}s'.format(read_count, t1 - t0))
  t0 = time.time()
  mio.sort_and_index_bam(args['<outbam>'])
  t1 = time.time()
  logger.debug('Sort and indexed perfect BAM in {:2.2f}s'.format(t1 - t0))


if __name__ == "__main__":
  cli()
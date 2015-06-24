"""In order to perform our alignment analysis we need to go through the BAM file and classify whether each read was
correctly aligned or not. The following back-of-the-envelope calculations indicate to us what size of data we are
expecting for a 30x illumina whole genome dataset:

# reads = 30 * 3e9 / 100 = 900e6 (900 million reads)

If about 5% if these are misaligned/unmapped, we have 45 million misaligned reads

For debugging purposes we store complete read data for only the incorrectly aligned reads. These are stored in a
database which offers efficient enough retrieval for debugging.

For comprehensive data analysis, where we want to compute alignment accuracy over different spots of the genome, but are
not interested in the read details - like sequence and base quality score - we only store the following information

chrom,
cpy,
start_pos  - 4 bytes  start of read
stop_pos   - 4 bytes  end of read    = start of read + read length
category   - 1 byte
---------------------
             9 bytes

The category flags are defined as follows

---------------------------------
| 7 | 6 | 5 | 4 | 3 | 2 | 1 | 0 |
---------------------------------
              |   |   |   |   |
              |   |   |   |   \------  1 => chrom was wrong
              |   |   |   \----------  1 => pos was wrong
              |   |   \--------------  1 => cigar was wrong
              |   \------------------  1 => unmapped
              \----------------------  1 => read from reference region (no variants)

To store this data for a 30x WG Illumina run will take 9 x 900e6 ~ 8.1 GB

This can  be loaded onto a modern machine. For this reason, we don't do anything very fancy and simply store the
categorized reads as a set of compressed numpy arrays using readily available mechanisms.

The categorized reads are stored in POS sorted order. The primary analysis we are are considering is to count the number
of correct and total reads in the vicinity of a series of positions (usually variants) or for the whole sequence without
regard to position (e.g. for reference reads).

The most efficient way to do this is to walk through the sorted list of reads side by side with the sorted list of
variants and count the relevant reads as we go.
"""
from mitty.version import __version__
__cmd__ = """perfectbam ({:s}): Categorize read alignments in a BAM file

Commandline::

  Usage:
    perfectbam <inbam> [--perfectbam=PBAM] [--debugdb=DB] [--catreads=CR] [--window=WN] [-x] [-v] [-p]

  Options:
    <inbam>             Input bam file name of reads
    --perfectbam=PBAM   Perfect BAM will be written to this file
    --debugdb=DB        Detailed data about misaligned reads will be written to this file
    --catreads=CR       All categorized reads will be saved to this file
    --window=WN         Size of tolerance window [default: 0]
    -x                  Use extended CIGAR ('X's and '='s) rather than traditional CIGAR (just 'M's)
    -v                  Dump detailed logger messages
    -p                  Show progress bar
""".format(__version__)

__param__ = """Given a bam file containing simulated reads aligned by a tool
  1. Produce a new bam that re-aligns all reads so that their alignment is perfect

  2. Produce a database file containing data about the misaligned and unmapped reads with the following tables::

  summary::

      chrom           - chromosome number
      total_reads     - total reads on this chrom
      incorrect_reads - total incorrect reads on this chrom
      unmapped_reads  - unmapped reads
      seq_len         - len of this sequence
      seq_id          - full seq name as found in BAM header

  reads::

       qname              -  qname of the read
       error_type         -  type of error 3 bit number  bit 0=chrom, 1=pos, 2=cigar, 3=unmapped
       correct_chrom      -  correct chromosome number of read
       correct_pos        -  correct position of read
       correct_cigar      -  correct cigar of read
       aligned_chrom      -  actual aligned chromosome
       aligned_pos        -  actual aligned pos
       aligned_cigar      -  actual aligned cigar
       mapping_qual       -  mapping quality
       mate_is_unmapped   -  is mate unmapped
       seq                -  actual sequence string

  3. Produce a compressed numpy array file with details of all the reads in the file

      chrom_cc   - 8 bit int: bit 0 -> chrom copy 0=0, 1=1  bit 1 onwards -> original chromosome of read
      pos        - 32 bit int: original position of read
      code       - reference_read_bit | error_type
                   bit 4 = 1 if read is from reference

      We could have split up the arrays into no_of_chroms x 2 arrays and dropped the chromosome field for space savings
      but this makes the code more complex since we no longer no how large to make the arrays.
      It also makes computations on the whole genome faster.
  """
__doc__ = __cmd__ + '\n\nDetails:\n\n' + __doc__

import os
import sqlite3 as sq
import time
import array

import numpy as np
import pysam
from collections import Counter
import docopt

import mitty.lib.io as mio  # For the bam sort and index function
from mitty.lib.reads import old_style_cigar
from mitty.lib import progress_bar

import logging
logger = logging.getLogger(__name__)


class CategorizedReads:
  """This class is a convenient interface to storing/analyzing read misalignment data.

  * In order to do our read alignment analysis we need access to the start and stop positions of each read


  * It wraps the error codes so we don't have to mess with/keep remembering bit fields.
  * Takes care of storing/reading the data from numpy arrays on file

  An example usage of this class is:

  raa = CategorizedReads(sequences=sequences, copies=copies)   #
  <in a loop, probably>
    raa.append(chrom, cpy, pos, code)
  raa.write(fname)

  You can load a previously saved analysis as:

  raa = CategorizedReads(fname='raa.npz')

  Internally, the data is stored as a set of C x c compressed numpy arrays where C is the number of chromosomes with
  c copies per chromosome.
  """
  def __init__(self, fname=None, seq_names=[], seq_lengths=[], copies=2):
    """

    :param fname:        If this is supplied an existing analysis is loaded from disk. All other params are ignored
    :param seq_names:    List of sequence ids
    :param seq_lengths:  Corresponding list of seq lengths
    :param copies:       How many copies per chromosome
    """
    if fname is not None:
      self.sequences = self.load_from_file(fname)
      self.copies = len(self.sequences[0][2])
      self.finalized = True
    else:
      assert len(seq_names) > 0
      assert len(seq_names) == len(seq_lengths)
      # Since there is no way to force a python array to be a specific width, we have to figure this out for ourselves
      dtypes = ['H', 'I', 'L']
      byte_size = [n for n, dt in enumerate(dtypes) if array.array(dt).itemsize == 4]
      if len(byte_size) == 0: raise ArithmeticError("Can't find an integer `array` type 4 bytes wide")
      dtype = dtypes[byte_size[0]]
      self.sequences = [[seq, l, [[array.array(dtype), array.array(dtype), array.array('B')] for _ in xrange(copies)]] for seq, l in zip(seq_names, seq_lengths)]
      self.copies = copies
      self.finalized = False

  def append(self, chrom, cpy, pos, read_len, code):
    """

    :param chrom:    correct chromosome number [1, 2, 3 ....]
    :param cpy:      correct copy number [0, 1, ...]
    :param pos:      correct pos
    :param read_len: length of read
    :param code:     code of read
    """
    a = self.sequences[chrom - 1][2][cpy]
    a[0].append(pos)
    a[1].append(pos + read_len)
    a[2].append(code)

  def write(self, fname):
    """Finalize and write data to numpy compressed file

    :param fname:
    """
    arrays_to_save = {
      'seq_names': np.array([s[0] for s in self.sequences]),
      'seq_lengths': np.array([s[1] for s in self.sequences]),
      'copies': len(self.sequences[0][2])
    }
    for chr_no, s in enumerate(self.sequences):
      for cpy, d in enumerate(s[2]):
        array_name = 'chrom_{:d}_copy_{:d}'.format(chr_no + 1, cpy)
        arrays_to_save[array_name + '_pos'] = np.frombuffer(d[0], dtype='uint32')
        arrays_to_save[array_name + '_stop'] = np.frombuffer(d[1], dtype='uint32')
        arrays_to_save[array_name + '_cat'] = np.frombuffer(d[2], dtype='B')
        srt_idx = np.argsort(arrays_to_save[array_name + '_pos'])
        #from IPython import embed; embed()
        arrays_to_save[array_name + '_pos'] = arrays_to_save[array_name + '_pos'][srt_idx]
        arrays_to_save[array_name + '_stop'] = arrays_to_save[array_name + '_stop'][srt_idx]
        arrays_to_save[array_name + '_cat'] = arrays_to_save[array_name + '_cat'][srt_idx]

    np.savez_compressed(fname, **arrays_to_save)
    self.finalized = True

  @staticmethod
  def load_from_file(fname):
    with np.load(fname) as fp:
      seq_names = fp['seq_names']
      seq_lengths = fp['seq_lengths']
      sequences = [
        [seq_name, seq_len,
          [
            [fp['chrom_{:d}_copy_{:d}_pos'.format(ch + 1, cpy)],
             fp['chrom_{:d}_copy_{:d}_stop'.format(ch + 1, cpy)],
             fp['chrom_{:d}_copy_{:d}_cat'.format(ch + 1, cpy)]] for cpy in range(fp['copies'])
          ]
        ]
        for ch, (seq_name, seq_len) in enumerate(zip(seq_names, seq_lengths))
      ]
    return sequences

  def get_data(self, chrom, cpy):
    assert 0 < chrom <= len(self.sequences), 'Enter a valid chromosome number'
    s = self.sequences[chrom - 1][2][cpy]
    return {'pos': s[0], 'stop': s[1], 'cat': s[2]}

  def __len__(self):
    return len(self.sequences)

  def chromosomes(self):
    return range(1, len(self) + 1)

  def __repr__(self):
    """Print out some useful information about ourselves"""
    rep_str = 'Read alignment analysis\n'
    rep_str += '{:d} chromosomes, {:d} copies each\n'.format(len(self.sequences), len(self.sequences[0][2]))
    for n in self.chromosomes():
      rep_str += 'Chrom {:d}: '.format(n)
      for c in range(self.copies):
         rep_str += 'copy {:d}: {:d} reads, '.format(c, self.get_data(n, c)['cat'].shape[0])
      rep_str += '\n'
    # TODO add more data
    return rep_str


def connect_to_db(db_name):
  """Convenience function so we don't need to import sqlite in other modules just for this"""
  conn = sq.connect(db_name)
  conn.row_factory = sq.Row
  return conn


def create_db(conn):
  """Create tables for the mis-aligned reads database

  :param conn: database connection

  read_serial = read_serial * 10 + 0 or 1 (for mate1 or mate2 of read) for paired reads
  read_serial = read_serial for un-paired reads
  """
  c = conn.cursor()
  c.execute('CREATE TABLE reads (read_serial INT, qname TEXT, error_type INT, '
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


def write_read_to_misaligned_read_table(conn, data_to_save):
  """Write mis-aligned/unmapped read data to the database

  :param conn: db connection
  :param data_to_save: [qname, error_type, chrom, pos, cigar, aligned_chrom, aligned_pos, aligned_cigar,
                        map_quality, mate_is_unmapped, sequence]
  """
  insert_clause = "INSERT INTO reads VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
  conn.execute(insert_clause, data_to_save)


def write_read_to_all_reads_table(conn, chrom_cc, pos, code):
  """Write mis-aligned/unmapped read data to the database

  :param conn: db connection
  :param chrom_cc:  correct chromosome of read & chrom copy bit
  :param pos: correct pos of read
  :param code: error code (same as error_type)
  """
  insert_clause = "INSERT INTO all_reads VALUES (?, ?, ?)"
  conn.execute(insert_clause, [chrom_cc, pos, code])


def commit_and_create_db_indexes(conn):
  """Need to commit the INSERTs and create search indexes for faster querying

  :param conn: db connection
  """
  conn.commit()
  conn.execute('CREATE INDEX idx_s ON reads (read_serial)')
  conn.execute('CREATE INDEX idx_f ON reads (correct_chrom, correct_pos)')
  conn.execute('CREATE INDEX idx_r ON reads (aligned_chrom, aligned_pos)')


def main(bam_in_fp, bam_out_fp=None, db_name=None, cr_fname=None, window=0, extended=False, progress_bar_func=None):
  """Main processing function that goes through the bam file, analyzing read alignment and writing out

  :param bam_in_fp:  Pointer to original BAM
  :param bam_out_fp: Pointer to perfect BAM being created
  :param db_name:    database name
  :param cr_fname:   name of file to store analysed reads in
  :param window:     Tolerance window for deciding if read is correctly aligned
  :param extended:   If True write out new style CIGARs (With '=' and 'X')
  :param progress_bar_func: Our famous progress bar, if we want it
  :return: number of reads processed
  """
  if db_name:
    conn = sq.connect(db_name)
    create_db(conn)
  else:
    conn = None

  total_read_count = float(bam_in_fp.mapped + bam_in_fp.unmapped)
  f0 = 0
  total_reads_cntr, incorrectly_aligned_reads_cntr, unmapped_reads_cntr = Counter(), Counter(), Counter()

  if progress_bar_func is not None: progress_bar_func('Processing BAM ', 0, 80)

  cat_reads = CategorizedReads(seq_names=bam_in_fp.references, seq_lengths=bam_in_fp.lengths) if cr_fname else None

  for n, read in enumerate(bam_in_fp):
    # qname = '{:d}|{:d}|{:d}|{:d}|{:s}|{:d}|{:s}'
    try:
      if read.is_paired:
        if read.is_read1:
          rs, chrom, cpy, ro, pos, cigar, _, _, _ = read.qname.split('|')
        else:
          rs, chrom, cpy, _, _, _, ro, pos, cigar = read.qname.split('|')
        read_serial = int(rs) * 10 + (not read.is_read1)
      else:
        rs, chrom, cpy, ro, pos, cigar = read.qname.split('|')[:6]  # For Wan-Ping :)
        read_serial = int(rs)
      ro, chrom, cpy, pos = int(ro), int(chrom), int(cpy), int(pos)
    except ValueError:
      logger.debug('Error processing qname: n={:d}, qname={:s}, chrom={:d}, pos={:d}'.format(n, read.qname, read.reference_id + 1, read.pos))
      continue

    # Do this before modifying the cigar
    error_type = 0b0 if 'X' in cigar or 'I' in cigar or 'D' in cigar or 'S' in cigar else 0b10000

    if not extended:
      cigar = old_style_cigar(cigar)

    total_reads_cntr[chrom] += 1

    if read.is_unmapped:
      error_type |= 0x8
    else:
      if read.reference_id != chrom - 1:
        error_type |= 1
      if not (-window <= read.pos - pos <= window):
        error_type |= 2
      if read.cigarstring != cigar:
        error_type |= 4
    if error_type & 0xff != 0:
      if error_type & 8:  # Unmapped read
        unmapped_reads_cntr[chrom] += 1
      else:
        incorrectly_aligned_reads_cntr[chrom] += 1
      if conn:
        write_read_to_misaligned_read_table(
          conn, [read_serial, read.qname, error_type, chrom, pos, cigar,
                 read.reference_id + 1, read.pos, read.cigarstring,
                 read.mapq, read.mate_is_unmapped, read.query_sequence])

    if cat_reads: cat_reads.append(chrom, cpy, pos, read.query_length, error_type)

    if bam_out_fp:
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

  if conn:
    write_summary_to_db(conn, total_reads_cntr, incorrectly_aligned_reads_cntr, unmapped_reads_cntr, bam_in_fp.header['SQ'])
    commit_and_create_db_indexes(conn)

  if cat_reads: cat_reads.write(cr_fname)

  return int(total_read_count)


def load_summary(conn):
  """Load information in the summary table

  :param conn: database connection object
  :return:
  """
  return [{k: row[k] for k in row.keys()} for row in conn.execute('SELECT * FROM summary ORDER BY chrom')]


def load_reads(conn, chrom=-1, start_pos=None, stop_pos=None, source=True, sample=None):
  """Load mis-aligned reads from the database based on the filters we give it

  :param conn: db connection object
  :param chrom: chromosome number [1,2...]. If -1, load all chroms and ignore start_pos and stop_pos values
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
  if chrom > 0:
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

  if args['--perfectbam'] is None and args['--debugdb'] is None and args['--catreads'] is None:
    print('No outputs specified. Easiest gig ever.')
    return

  bam_in_fp = pysam.AlignmentFile(args['<inbam>'], 'rb')
  pbam_fname = args['--perfectbam']
  bam_out_fp = pysam.AlignmentFile(pbam_fname, 'wb', template=bam_in_fp) if pbam_fname else None
  db_name = args['--debugdb']
  cr_fname = args['--catreads']
  if db_name:
    try:
      os.remove(db_name)
    except OSError:
      pass

  t0 = time.time()
  read_count = main(bam_in_fp=bam_in_fp, bam_out_fp=bam_out_fp, db_name=db_name, cr_fname=cr_fname,
                    window=int(args['--window']), extended=bool(args['-x']), progress_bar_func=progress_bar if args['-p'] else None)
  if bam_out_fp: bam_out_fp.close()
  t1 = time.time()
  logger.debug('Analyzed {:d} reads in BAM in {:2.2f}s'.format(read_count, t1 - t0))

  if pbam_fname:
    t0 = time.time()
    mio.sort_and_index_bam(args['--perfectbam'])
    t1 = time.time()
    logger.debug('Sort and indexed perfect BAM in {:2.2f}s'.format(t1 - t0))


if __name__ == "__main__":
  cli()
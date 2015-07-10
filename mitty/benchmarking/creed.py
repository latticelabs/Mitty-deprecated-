"""Functions to categorize reads for further analysis will be refactored into this file.




"""
__file_formats__ = """The bracket file should be an HDF5 file. Each bracket array is a structured numpy array with
dtype=[('pos', 'i4'), ('stop', 'i4')]. The outputs are structured arrays ('buckets')
dtype=[('correct', 'i4'), ('total', 'i4')].
There is one bucket array for each bracket array and each bucket array is the same length as the corresponding bucket
array."""

import sqlite3 as sq
from itertools import izip
import array

import h5py
import numpy as np
from pysam import AlignedSegment as pas

from mitty.lib.reads import old_style_cigar

import logging
logger = logging.getLogger(__name__)


class CategorizedReads:
  """A convenient interface to storing/retrieving read misalignment data. Stores the start, end and category of each
  read in the bam file

  cat_reads = CategorizedReads(fname='catreads.hdf5', seq_names=bam_in_fp.references, seq_lengths=bam_in_fp.lengths)
  <in a loop, probably>
    raa.append(chrom, cpy, pos, code)
  raa.write(fname)

  You can load a previously saved analysis as:

  raa = CategorizedReads(fname='raa.npz')

  Internally, the data is stored as a set of C x c compressed numpy arrays where C is the number of chromosomes with
  c copies per chromosome.
  """
  def __init__(self, fname=None, seq_names=[], seq_lengths=[], copies=None):
    """

    :param fname:        File name to save/load data. If None, create a temporary file in memory
    :param seq_names:    List of sequence ids              # If this and the following arguments are given a new
    :param seq_lengths:  Corresponding list of seq lengths # file is started. Otherwise an existing file is loaded
    :param copies:       How many copies per chromosome
    """
    if fname is None:
      self.fp = h5py.File(name=hex(id(self)), driver='core', backing_store=False)
    else:
      self.fp = h5py.File(name=fname)

    if len(self.get_sequence_info()) == 0:  # We are not loading a file but creating it
      if len(seq_names) == 0 or len(seq_lengths) == 0 or copies is None:
        raise RuntimeError('Creating a new CategorizedReads object requires genome metadata and chrom copy count')
      assert len(seq_names) == len(seq_lengths)
      self.fp.attrs['sequence_names'] = seq_names
      self.fp.attrs['sequence_lens'] = seq_lengths
      self.fp.attrs['copies'] = copies
      i4 = self.get_i4type()
      self._py_arrays = [
        [[array.array(i4), array.array(i4), array.array('B')] for _ in range(copies)]
        for _ in range(len(seq_names))
      ]
    self.copies = self.fp.attrs['copies']

  @staticmethod
  def get_i4type():
    # Since there is no way to force a python array to be a specific width, we have to figure this out for ourselves
    dtypes = ['H', 'I', 'L']
    byte_size = [n for n, dt in enumerate(dtypes) if array.array(dt).itemsize == 4]
    if len(byte_size) == 0: raise ArithmeticError("Can't find an integer `array` type 4 bytes wide")
    return dtypes[byte_size[0]]

  def get_sequence_info(self):
    """
    :return: [(seq_id, seq_len) ...]
    """
    return [(n, l) for n, l in zip(self.fp.attrs.get('sequence_names', []), self.fp.attrs.get('sequence_lens', []))]

  def append(self, chrom, cpy, pos, read_len, code):
    """

    :param chrom:    correct chromosome number [1, 2, 3 ....]
    :param cpy:      correct copy number [0, 1, ...]
    :param pos:      correct pos
    :param read_len: length of read
    :param code:     code of read
    """
    try:
      a = self._py_arrays[chrom - 1][cpy]
      a[0].append(pos)
      a[1].append(pos + read_len)
      a[2].append(code)
    except AttributeError:
      raise RuntimeError("Can't append data after finalizing")

  def finalize(self):
    """Convert py arrays to numpy structured arrays in HDF5 format. Remove _py_arrays so we can no longer append"""
    dtype = [('pos', 'i4'), ('stop', 'i4'), ('cat', 'B')]
    for n in range(len(self.fp.attrs['sequence_names'])):
      for cpy in range(self.copies):
        temp_array = np.empty(shape=(len(self._py_arrays[n][cpy][0]),), dtype=dtype)
        for dt in [0, 1, 2]:
          temp_array[dtype[dt][0]] = self._py_arrays[n][cpy][dt]
        dset = self.fp.create_dataset(name='/chrom_{:d}/copy_{:d}'.format(n + 1, cpy),
                                      shape=(len(self._py_arrays[n][cpy][0]),),
                                      data=np.sort(temp_array, order=['pos']), dtype=dtype,
                                      chunks=True, compression='gzip')
        dset.attrs['sequence_name'] = self.fp.attrs['sequence_names'][n]
    del self._py_arrays

  def get_data(self, chrom, cpy):
    assert 0 < chrom <= len(self.fp.attrs['sequence_names']), 'Enter a valid chromosome number'
    return self.fp['/chrom_{:d}/copy_{:d}'.format(chrom, cpy)][:]

  def __len__(self):
    return len(self.fp.attrs['sequence_names'])

  def chromosomes(self):
    return range(1, len(self) + 1)

  def __repr__(self):
    """Print out some useful information about ourselves"""
    rep_str = 'Read alignment analysis\n'
    rep_str += '{:d} chromosomes, {:d} copies each\n'.format(len(self), self.copies)
    for n in self.chromosomes():
      rd_cnts = ['{:d}'.format(self.fp['/chrom_{:d}/copy_{:d}'.format(n, c)].shape[0]) for c in range(self.copies)]
      rep_str += 'Chrom {:d} ['.format(n) + '|'.join(rd_cnts) + ']'
      rep_str += '\n'
    # TODO add more data
    return rep_str


class ReadDebugDB:
  def __init__(self, db_name):
    """Create a new database with this name"""
    self.db_name = db_name
    self.conn = sq.connect(self.db_name)
    self.conn.row_factory = sq.Row
    if not self.db_exists():
      self.create_db()

  def db_exists(self):
    if len([r for r in self.conn.execute('PRAGMA table_info(reads)')]):
      logger.debug('Database exists')
      return True
    else:
      return False

  def create_db(self):
    """Create tables for the mis-aligned reads database

    read_serial = read_serial * 10 + 0 or 1 (for mate1 or mate2 of read) for paired reads
    read_serial = read_serial for un-paired reads
    """
    c = self.conn.cursor()
    c.execute('CREATE TABLE reads (read_serial INT, qname TEXT, error_type INT, '
              'correct_chrom INT, correct_pos INT, correct_cigar TEXT, '
              'aligned_chrom INT, aligned_pos INT, aligned_cigar TEXT,'
              'mapping_qual INT, mate_is_unmapped BOOL, seq TEXT)')
    c.execute('CREATE TABLE summary (chrom INT, total_reads INT, incorrect_reads INT, unmapped_reads INT, seq_len INT, seq_id TEXT)')
    self.conn.commit()

  def write_summary_to_db(self, total_reads_cntr, incorrectly_aligned_reads_cntr, unmapped_reads_cntr, bam_seq_header):
    """Write the alignment analysis summary for the BAM file into the database

    :param total_reads_cntr: Counter by chromosome
    :param incorrectly_aligned_reads_cntr: Counter by chromosome
    :param unmapped_reads_cntr: Counter by chromosome
    :param bam_seq_header: bam_in_fp.header['SQ']
    :return:
    """
    insert_clause = "INSERT INTO summary VALUES (?, ?, ?, ?, ?, ?)"
    to_insert = [(n + 1, total_reads_cntr[n + 1], incorrectly_aligned_reads_cntr[n + 1], unmapped_reads_cntr[n + 1], seq['LN'], seq['SN']) for n, seq in enumerate(bam_seq_header)]
    self.conn.executemany(insert_clause, to_insert)

  def write_read_to_misaligned_read_table(self, data_to_save):
    """Write mis-aligned/unmapped read data to the database

    :param data_to_save: [qname, error_type, chrom, pos, cigar, aligned_chrom, aligned_pos, aligned_cigar,
                          map_quality, mate_is_unmapped, sequence]
    """
    insert_clause = "INSERT INTO reads VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
    self.conn.execute(insert_clause, data_to_save)

  def commit_and_create_db_indexes(self):
    """Need to commit the INSERTs and create search indexes for faster querying"""
    conn = self.conn
    conn.commit()
    conn.execute('CREATE INDEX idx_s ON reads (read_serial)')
    conn.execute('CREATE INDEX idx_f ON reads (correct_chrom, correct_pos)')
    conn.execute('CREATE INDEX idx_r ON reads (aligned_chrom, aligned_pos)')

  def load_summary(self):
    """Load information in the summary table

    :return:
    """
    return [{k: row[k] for k in row.keys()} for row in self.conn.execute('SELECT * FROM summary ORDER BY chrom')]

  def load_reads(self, chrom=-1, start_pos=None, stop_pos=None, source=True, sample=None):
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
    return [r for r in self.conn.execute(query)]


def analyze_read(read, window=100, extended=False):
  """Given a read process the qname and read properties to determine what kind of alignment errors were made on it and

  :param read: a psyam AlignedSegment object
  :returns read_serial, chrom, cpy, ro, pos, cigar, read_category


  read_serial = read_serial * 10 + 0 or 1 (for mate1 or mate2 of read) for paired reads
  read_serial = read_serial for un-paired reads

  The category flags are defined as follows

  ---------------------------------
  | 7 | 6 | 5 | 4 | 3 | 2 | 1 | 0 |
  ---------------------------------
            |   |   |   |   |   |
            |   |   |   |   |   \------  1 => chrom was wrong
            |   |   |   |   \----------  1 => pos was wrong
            |   |   |   \--------------  1 => cigar was wrong
            |   |   \------------------  1 => unmapped
            |   \----------------------  1 => read from reference region (no variants)
            \--------------------------  1 => mate is from reference region
  """
  early_exit_value = [None] * 7

  # Not counted
  if read.is_secondary:
    return early_exit_value

  # We should never actually fail this, unless a tool messes badly with the qname
  try:
    if read.is_paired:
      if read.is_read1:
        rs, chrom, cpy, ro, pos, cigar, _, _, cigar1 = read.qname.split('|')  # We decode mate cigar to determine if mate is a reference read
      else:
        rs, chrom, cpy, _, _, cigar1, ro, pos, cigar = read.qname.split('|')
      read_serial = int(rs) * 10 + (not read.is_read1)
    else:
      rs, chrom, cpy, ro, pos, cigar = read.qname.split('|')[:6]  # For Wan-Ping :)
      read_serial = int(rs)
    ro, chrom, cpy, pos = int(ro), int(chrom), int(cpy), int(pos)
  except ValueError:
    logger.debug('Error processing qname: qname={:s}, chrom={:d}, pos={:d}'.format(read.qname, read.reference_id + 1, read.pos))
    return early_exit_value

  # Do this before modifying the cigar
  read_category = 0b000000
  if not ('X' in cigar or 'I' in cigar or 'D' in cigar or 'S' in cigar):
    read_category |= 0b010000
  if read.is_paired:
    if not ('X' in cigar1 or 'I' in cigar1 or 'D' in cigar1 or 'S' in cigar1):
      read_category |= 0b100000

  if not extended:
    cigar = old_style_cigar(cigar)

  if read.is_unmapped:
    read_category |= 0b1000
  else:
    if read.reference_id != chrom - 1:
      read_category |= 0b0001
    # if not (-window <= read.pos - pos <= window):
    #   read_category |= 0b0010
    perfect_read = pas()
    perfect_read.cigarstring = cigar
    perfect_read.pos = pos
    if not (perfect_read.positions[0] <= read.pos <= perfect_read.positions[-1]):
      read_category |= 0b0010
    if read.cigarstring != cigar:
      read_category |= 0b0100

  return read_serial, chrom, cpy, ro, pos, cigar, read_category


def bucket_list(r_pos, r_stop, r_code, v_pos, v_stop):
  """First run - ignoring CIGAR errors

  :param r_pos:  array of read 'pos'  -> should be sorted
  :param r_stop: array of read 'stop'
  :param r_code: array of read category code
  :param v_pos: start of bucket  -> should be sorted
  :param v_stop: end of bucket
  :return:
  """
  # correct, total
  v_win_start = 0  # Which variants do we check
  v_win_stop = 0
  v_count = v_pos.shape[0]
  v_count_1 = v_count - 1
  v_read_counts = np.zeros(v_count, dtype=[('correct', 'uint32'), ('total', 'uint32')])
  v_correct = v_read_counts['correct']
  v_total = v_read_counts['total']
  ref_one_correct, ref_one_total = 0, 0
  ref_both_correct, ref_both_total = 0, 0  # Both read and mate are from reference regions

  for r0, r1, c in izip(r_pos, r_stop, r_code):
    # If this is a reference read we can bucket it now and move on
    if c & 0b10000:  # reference read
      if c & 0b100000:  # mate is reference
        ref_both_total += 1
        if not (c & 0b1011):
          ref_both_correct += 1
      else:  # mate is not reference
        ref_one_total += 1
        if not (c & 0b1011):
          ref_one_correct += 1
      continue

    if v_count == 0: continue  # No variants, don't bother

    # This is a variant read and we have variants

    # Advance our variant window as needed
    while v_stop[v_win_start] < r0:
      if v_win_start < v_count_1:
        v_win_start += 1
      else:
        break

    if v_win_stop < v_count_1:  # Only need to advance if there is something left
      while v_pos[v_win_stop + 1] < r1:
        v_win_stop += 1
        if v_win_stop == v_count_1:
          break

    for n in range(v_win_start, v_win_stop + 1):
      if r0 <= v_stop[n] and r1 >= v_pos[n]:
        v_total[n] += 1
        if not (c & 0b1011):  # no chrom, pos or unmapped errors
          v_correct[n] += 1

  return np.array([[ref_one_correct, ref_one_total], [ref_both_correct, ref_both_total]], dtype='i4'), v_read_counts


def categorize_read_counts_by_indel_length(variations, v_read_counts, cat_counts=None, max_indel=100):
  assert max_indel > 0
  assert len(variations) == len(v_read_counts)
  if cat_counts is None:
    cat_counts = np.zeros(2 * max_indel + 1, dtype=[('x', 'int32'), ('correct', 'uint32'), ('total', 'uint32')])
    cat_counts['x'] = range(-max_indel, max_indel + 1)  # The range of indel lengths we are assuming.
  correct = cat_counts['correct']
  total = cat_counts['total']
  for v, c in izip(variations, v_read_counts):
    d = len(v['alt']) - len(v['ref'])
    if abs(d) > max_indel: continue
    total[d + max_indel] += c['total']
    correct[d + max_indel] += c['correct']
  return cat_counts


def find_nearest_variant(v1, v2):
  """For every variant in v1 find the nearest variant in v2 and note its distance and length. length = alt - ref such
  that SNP = 0, insertion > 0 and deletion < 0

  :param v1: list 1
  :param v2: list 2
  :return: numpy structured array the same length as v1 with fields 'dist' and 'length'
  """
  nearest_v2 = np.empty(v1.shape[0], dtype=[('dist', 'uint32'), ('length', 'int32')])
  nearest_dist, nearest_len = nearest_v2['dist'], nearest_v2['length']
  nearest_dist[:], nearest_len[:] = np.Inf, np.nan

  v2_size = v2.shape[0]
  v2_size_1 = v2_size - 1
  if v2_size == 0: return nearest_v2
  v2_idx = 0
  for n, v in enumerate(v1):
    d_min = abs(v2[v2_idx]['pos'] - v['pos'])
    while 1:
      if v2_idx < v2_size_1:
        d = abs(v2[v2_idx + 1]['pos'] - v['pos'])
        if d <= d_min:
          d_min = d
          v2_idx += 1
      else:
        break
    nearest_dist[n] = d_min
    nearest_len[n] = v2[n]['alt'] - v2[n]['ref']
  return nearest_v2


def categorize_read_counts_by_indel_length_and_nearest_variant(v1, v_read_counts, nearest_v2,
                                                               cat_counts=None,
                                                               max_v1_indel=100,
                                                               max_v2_indel=100,
                                                               max_dist=300):
  """v2 indels that are beyond the max range are piled into the first or last bins."""
  assert max_v1_indel > 0
  assert max_v2_indel > 0
  assert len(v1) == len(v_read_counts)
  assert len(v1) == len(nearest_v2)
  if cat_counts is None:
    cat_counts = {
      'counts': np.zeros((2 * max_v1_indel + 1, 2 * max_v2_indel + 1, max_dist + 1), dtype=[('correct', 'uint32'), ('total', 'uint32')]),
      'indel_v1_size': range(-max_v1_indel, max_v1_indel + 2),
      'indel_v2_size': range(-max_v1_indel, max_v1_indel + 1),
      'dist': range(max_dist + 1)
    }
  correct = cat_counts['counts']['correct']
  total = cat_counts['counts']['total']
  for n in xrange(v1.shape[0]):
    d_v1 = len(v1[n]['alt']) - len(v1[n]['ref'])
    if abs(d_v1) > max_v1_indel: continue
    idx1 = max_v1_indel + d_v1
    idx2 = max(0, min(2 * max_v2_indel, nearest_v2[n]['length'] + max_v2_indel))
    idx3 = min(max_dist + 1, nearest_v2[n]['dist'])
    total[idx1, idx2, idx3] += v_read_counts['total']
    correct[idx1, idx2, idx3] += v_read_counts['correct']
  return cat_counts
"""Library to handle database interactions.

The master list of variants are stored in M tables named master_chrom_1, master_chrom_2, ... master_chrom_M

Columns::

  rowid  - corresponds to index
  pos    -
  stop   -
  ref    -
  alt    -
  p      - probability of variant (to make the sfs)


Sample data is stored in M tables named sample_chrom_1, ....

Columns::

  rowid  - sample id
  gen    - generation
  serial - serial within the generation
  data   - binary blob

          Chromosome data is stored as a sequence of 4N bytes
          30 bits represent the index and 2 bits represent the genotype
          This is stored as a blob in the vlist field

"""
import sqlite3 as sq
import struct

import numpy

import mitty.lib.variants as vr


def sample_table_name(sample_name, chrom_name):
  return 'sample_{:s}_chrom_{:s}'.format(sample_name, chrom_name)


def connect(db_name='population.sqlite3'):
  """Connect to the database
  :param db_name: The database name
  :returns conn: connection object"""
  conn = sq.connect(db_name)
  conn.text_factory = str  # Otherwise our ref and alts will be unicode, bleh!
  return conn


def save_master_list(conn, chrom, ml):
  """
  :param conn: connection object
  :param chrom: chromosome number
  :param ml: master list of variants
  :return: nothing

  THIS WIPES ANY PREVIOUS DATA AND WRITES THE LIST AFRESH
  """
  assert type(chrom) == int, 'Chromosome must be a number'
  assert chrom > 0, 'Chromosome numbering starts at 1'
  assert ml.sorted, 'Master list has not been sorted. Please check your program'
  assert len(ml) <= 1073741823, 'Master list has more than 2^30-1 variants.'  # I want whoever gets here to mail me: kaushik.ghose@sbgenomics.com

  conn.execute("DROP TABLE IF EXISTS master_chrom_{:d}".format(chrom))
  conn.execute("CREATE TABLE master_chrom_{:d} (pos INTEGER, stop INTEGER, ref TEXT, alt TEXT, p FLOAT)".format(chrom))

  insert_clause = "INSERT INTO master_chrom_{:d} (pos, stop, ref, alt, p) VALUES (?, ?, ?, ?, ?)".format(chrom)
  conn.executemany(insert_clause, ml.variants.tolist())
  conn.commit()


def load_master_list(conn, chrom):
  """
  :param conn: connection object
  :param chrom: chromosome number
  :returns ml: master list of variants
  """
  assert type(chrom) == int, 'Chromosome must be a number'
  assert chrom > 0, 'Chromosome numbering starts at 1'

  # # Surely numpy has a better way of doing this, but fromiter does not work
  # c = []
  # for col in ['pos', 'stop', 'ref', 'alt', 'p']:
  #   c += [[r for r in conn.execute("SELECT {:s} FROM master_chrom_{:d}".format(col, chrom))]]
  #ml = vr.VariantList(*c)

  ml = vr.VariantList()

  try:
    rows = [c for c in conn.execute("SELECT * FROM master_chrom_{:d}".format(chrom))]
    dtype = [('pos', 'i4'), ('stop', 'i4'), ('ref', 'object'), ('alt', 'object'), ('p', 'f2')]
    ml.variants = numpy.core.records.fromrecords(rows, dtype=dtype)
  except sq.OperationalError:
    pass  # TODO: log some kind of warning? Or assume we just don' have any variants for that chromosome?
  ml.sorted = True  # We assume that this had been sorted etc. before saving
  return ml


def save_sample(conn, gen, serial, chrom, variants):
  """Save the sample chromosomes into the database

  :param conn: connection object
  :param gen: generation
  :param serial: serial number within the generation
  :param chrom: chromosome number
  :param variants: list of tuples as returned by generate_chromosome
  :return: rowid corresponding to this sample

  Chromosome data is stored as a sequence of 4N bytes
  30 bits represent the index and 2 bits represent the genotype
  This is stored as a blob in the vlist field
  """
  c = conn.cursor()
  c.execute("CREATE TABLE IF NOT EXISTS sample_chrom_{:d} (gen INTEGER, serial INTEGER, vlist BLOB)".format(chrom))
  c.execute("INSERT INTO sample_chrom_{:d} (gen, serial, vlist) VALUES (?, ?, ?)".format(chrom), (gen, serial,
            sq.Binary(struct.pack('{:d}I'.format(len(variants)), *[v[0] << 2 | v[1] for v in variants]))))
  conn.commit()
  return c.lastrowid


def load_sample(conn, gen, serial, chrom):
  """Load the sample in the database

  :param conn: connection object
  :param gen: generation
  :param serial: serial number within the generation
  :param chrom: chromosome number
  :return: variants: list of tuples same as returned by generate_chromosome
  """
  c = conn.cursor()
  try:
    c.execute("SELECT vlist FROM sample_chrom_{:d} WHERE gen==? AND serial==?".format(chrom), (gen, serial))
    row = next(c, None)
  except sq.OperationalError:
    row = None  # We assume that no such table means, simply, no variants in that chromosome
  return [(b >> 2, b & 0x3) for b in struct.unpack('{:d}I'.format(len(row[0]) / 4), row[0])] if row is not None else []


def save_chromosome_metadata(conn, chrom, seq_id, seq_len, seq_md5):
  """Save chromosome sequence metadata in db

  :param conn:    connection object
  :param chrom:   chromosome number
  :param seq_id:  sequence id as found in fasta
  :param seq_len: sequence length
  :param seq_md5:md5 hash of sequence string
  """
  c = conn.cursor()
  c.execute("CREATE TABLE IF NOT EXISTS chrom_metadata (chrom INTEGER, seq_id TEXT, seq_len INTEGER, seq_md5 TEXT)")
  c.execute("INSERT INTO chrom_metadata (chrom, seq_id, seq_len, seq_md5) VALUES (?, ?, ?, ?)", (chrom, seq_id, seq_len, seq_md5))
  conn.commit()


def load_chromosome_metadata(conn, chrom):
  """Load chromosome metadata

  :param conn:   connection object
  :param chrom:  chromosome number
  :returns chrom, seq_id, seq_len, seq_md5
  """
  c = conn.cursor()
  c.execute("SELECT chrom, seq_id, seq_len, seq_md5 FROM chrom_metadata WHERE chrom=?", (chrom,))
  row = next(c, None)
  return row


def chromosomes_in_db(conn):
  # c = conn.execute("SELECT name FROM sqlite_master WHERE TYPE='table' AND name LIKE 'master_chrom_%'")
  # return [int(row[0].replace('master_chrom_','')) for row in c]
  c = conn.execute("SELECT chrom, seq_id, seq_len, seq_md5 FROM chrom_metadata ORDER BY rowid ASC")
  return [row for row in c]


def variants_in_master_list(conn, chrom):
  c = conn.execute("SELECT COUNT(rowid) FROM master_chrom_{:d}".format(chrom))
  return c.next()[0]


def samples_in_db(conn):
  c = conn.execute("SELECT name FROM sqlite_master WHERE TYPE='table' AND name LIKE 'sample_chrom_%'")
  table = c.next()[0]
  c = conn.execute("SELECT COUNT(rowid) FROM {:s}".format(table))
  return int(c.next()[0])

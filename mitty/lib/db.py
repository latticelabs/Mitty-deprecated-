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

  dtype = [('pos', 'i4'), ('stop', 'i4'), ('ref', 'object'), ('alt', 'object'), ('p', 'f2')]
  rows = [c for c in conn.execute("SELECT * FROM master_chrom_{:d}".format(chrom))]
  ml = vr.VariantList()
  ml.variants = numpy.core.records.fromrecords(rows, dtype=dtype)
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
  c.execute("SELECT vlist FROM sample_chrom_{:d} WHERE gen==? AND serial==?".format(chrom), (gen, serial))
  row = next(c, None)
  return [(b >> 2, b & 0x3) for b in struct.unpack('{:d}I'.format(len(row[0]) / 4), row[0])] if row is not None else []
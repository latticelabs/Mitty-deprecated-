"""Library to handle database interactions.

tables in database (X is the chromosome number):

mitty.lib.variation.Variant -> var_X
mitty.lib.variation.Chromosome -> samp_I_chrom_X
mitty.lib.variation.Sample -> samples

var_X
rowid, pos, stop, ref, alt  rowid corresponds to index

sample_I_chrom_X
rowid, index, gt  rowid is in order of variant. index corresponds to var_X

samples
rowid, generation, serial, fitness, p1, p2    p1, p2 are the serials of the parents. We now the generation is one less
"""
import sqlite3 as sq
import mitty.lib.variation as vr


# def var_factory(cursor, row):
#   """For speed, assumes we load the columns in the order they are in the db."""
#   v = vr.new_variant(row[1], row[2], row[3], row[4])
#   v.index = row[0]
#   return v


def sample_table_name(sample_name, chrom_name):
  return 'sample_{:s}_chrom_{:s}'.format(sample_name, chrom_name)


def connect(db_name='population.sqlite3'):
  """Connect to the database
  :param db_name: The database name
  :returns conn: connection object"""
  conn = sq.connect(db_name)
  #conn.text_factory = str
  #conn.row_factory = var_factory
  return conn


def save_variant_master_list(conn, ml):
  """
  :param conn: connection object
  :param ml: master list of variants
  :return: nothing

  For now, we assume this is always a fresh database and we are writing all the variants in one go
  """
  for chrom in ml.keys():
    conn.execute("CREATE TABLE chrom_{:d} (pos INTEGER, stop INTEGER, ref TEXT, alt TEXT)".format(chrom))
    insert_clause = "INSERT INTO chrom_{:d} (rowid, pos, stop, ref, alt) VALUES (?, ?, ?, ?, ?)".format(chrom)
    conn.executemany(insert_clause, (v.as_tuple() for v in ml[chrom].values()))
  conn.commit()


def erase_sample(sample_name, g, conn):
  """Drop tables associated with this sample
  :param sample_name: name of the sample
  :param g: Genome object. Only needed for the chromosome ids
  :param conn: The connection object.
  """
  for chrom_name in g.keys():
    conn.execute("DROP TABLE IF EXISTS {:s}".format(sample_table_name(sample_name, chrom_name)))


def save_sample(sample_name, g, conn):
  """Save the sample in the database
  :param sample_name: name of the sample
  :param g: Genome object
  :param conn: The connection object.

  Note: this erases previous data for this sample
  """
  erase_sample(sample_name, g, conn)
  for k, chrom in g.iteritems():
    table_name = sample_table_name(sample_name, k)
    conn.execute("CREATE TABLE {:s} (pos INTEGER PRIMARY KEY, stop INTEGER, gt INTEGER, ref TEXT, alt TEXT)".format(table_name))
    insert_clause = "INSERT INTO {:s}(pos, stop, gt, ref, alt) VALUES (?, ?, ?, ?, ?)".format(table_name)
    conn.executemany(insert_clause, (c.as_tuple() for c in chrom))
  conn.commit()


def load_sample(sample_name, g, conn):
  """Load the sample in the database
  :param sample_name: name of the sample
  :param g: Genome object. Only needed for the chromosome ids. Any actual data in this object will be erased
  :param conn: The connection object.
  :returns g implicitly (changes g in place)
  """
  for k, chrom in g.iteritems():
    table_name = sample_table_name(sample_name, k)
    g[k] = [r for r in conn.execute("SELECT * FROM {:s}".format(table_name))]


# def initialize_population_table(conn):
#   """Create a new table to store population data.
#   DESTROYS ANY EXISTING POPULATION TABLE
#   :param conn: The connection object
#   """
#   conn.execute("DROP TABLE IF EXISTS population")
#   conn.execute("CREATE TABLE population (gen INTEGER, serial INTEGER, p1 INTEGER, p2 INTEGER, "
#                "FOREIGN KEY(p1) REFERENCES population(rowid), FOREIGN KEY(p2) REFERENCES population(rowid))")
#
#
# def save_subpopulation(conn):


"""Library to handle database interactions"""
import sqlite3 as sq

# Types of zygosity
ABSENT = 0
HOMOZYGOUS = 3
HET_01 = 1
HET_10 = 2


def sample_table_name(sample_name, chrom_name):
  return 'sample_{:s}_chrom_{:s}'.format(sample_name, chrom_name)


def connect(db_name='population.sqlite3'):
  """Connect to the database
  :param db_name: The database name
  :returns conn: connection object"""
  conn = sq.connect(db_name)
  conn.text_factory = str
  return conn


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
    conn.executemany(insert_clause, chrom)
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


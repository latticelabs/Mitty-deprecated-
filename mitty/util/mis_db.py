"""Useful functions to construct queries on the misalignment database(s) produced by perfectbam."""


def attach_databases(conn, db_name_list):
  """Attach the databases to the conn instance

  :param conn:           a db connection (dummy/empty ones, like ':memory:' work nicely)
  :param db_name_list:   list of database names
  """
  for n, db_name in enumerate(db_name_list):
    conn.execute('ATTACH DATABASE ? AS db{:d}'.format(n), (db_name,))


def from_clause(db1, db2=None, operation='diff', chrom=None, variants_only=False):
  """Create a 'from' clause from the parameters

  :param db1:        row [0, 1, 2 ...]
  :param db2:        col [0, 1, 2 ...]. Omit to do selects from db1 only
  :param operation:  'diff' or 'intersection' ('d' or 'i'). Ignored if db2 is none
  :param chrom:      chromosome number. Omit to do all chromosomes
  :param variants_only: Set to True to only pick reads with variants in them
  :return: a from clause
  """
  cc = ' AND t1.correct_chrom={:d}'.format(chrom) if chrom else ''
  vv = " AND (t1.correct_cigar LIKE '%I%' OR t1.correct_cigar LIKE '%S%' " \
       "OR t1.correct_cigar LIKE '%D%' OR t1.correct_cigar LIKE '%X%')" if variants_only else ''
  if db2 is not None:
    if operation.startswith('d'):
      # http://blog.montmere.com/2010/12/08/the-anti-join-all-values-from-table1-where-not-in-table2/
      clause = 'from db{db1}.reads as t1 left join db{db2}.reads as t2' \
               ' on t1.read_serial = t2.read_serial where t2.read_serial is NULL {cc} {vv};'
    elif operation.startswith('i'):
      clause = 'from db{db1}.reads as t1, db{db2}.reads as t2' \
               ' WHERE t1.read_serial = t2.read_serial {cc} {vv};'
    else:
      raise ValueError, 'No such option' + operation
  else:
      clause = 'from db{db1}.reads as t1 ' + 'WHERE' if cc or vv else '' + ' {cc} {vv};'

  return clause.format(db1=db1, db2=db2, cc=cc, vv=vv)


def get_bench_mat(conn, db_count, chrom=None, variants_only=False):
  """For all the valid combinations of row and col, fill out the benchmarking matrix.

  :param conn:       database connection with databases attached
  :param db_count:   number of attached databases
  :param chrom:      chromosome number. Omit to do all chromosomes
  :param variants_only: Set to True to only pick reads with variants in them
  :returns dictionary representing cells in benchmark matrix
  """
  bm_mat = {'metadata': {'chrom': chrom, 'dbs': {r[1]: r[2] for r in conn.execute('PRAGMA database_list') if r[1].startswith('db')}}}
  for r in range(db_count):
    for c in range(db_count):
      if c != r:
        bm_mat[(r, c, 'd')] = count_mat_cell(conn, r, c, operation='diff', chrom=chrom, variants_only=variants_only)
      if c >= r:
        bm_mat[(r, c, 'i')] = count_mat_cell(conn, r, c, operation='i', chrom=chrom, variants_only=variants_only)
  return bm_mat


def count_mat_cell(conn, db1, db2, operation='diff', chrom=None, variants_only=False):
  """Get read count for given cell and operation

  :param conn:       connection with databases attached
  :param db1:        row [0, 1, 2 ...]
  :param db2:        col [0, 1, 2 ...]
  :param operation:  'diff' or 'intersection' ('d' or 'i')
  :param chrom:      chromosome number. Omit to do all chromosomes
  :param variants_only: Set to True to only pick reads with variants in them
  :returns integer representing read count for this cell
  """
  c = conn.execute('SELECT COUNT(*) ' + from_clause(db1, db2, operation, chrom, variants_only))
  return c.next()[0]


def get_reads_mat_cell(conn, db1, db2=None, operation='diff', chrom=None, variants_only=False):
  """Counterpart of count_mat_cell that returns us the actual reads in this 'cell'

  :param conn:       connection with databases attached
  :param db1:        row [0, 1, 2 ...]
  :param db2:        col [0, 1, 2 ...]
  :param operation:  'diff' or 'intersection' ('d' or 'i')
  :param chrom:      chromosome number. Omit to do all chromosomes
  :param variants_only: Set to True to only pick reads with variants in them
  :returns Generator from database.
  """
  return conn.execute('SELECT * ' + from_clause(db1, db2, operation, chrom, variants_only))
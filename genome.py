"""This module defines a whole genome file (.wg) that allows us to store multiple sequences in an organized fashion. It
also defines an interface to read in chosen sequences from the file. The file is meant to represent an individual and
carries multiple chromosomes and multiple copies of each chromosome. The rationale for using this (rather than a pickle
or some other format) is that it is a simple file format that allows us to efficiently access any part of the genome.
We could have used HDF5 but it would have been standards compliant at the expense of unnecessary complexity.

Usage:
genome.py --wg=WG
genome.py explain

Options:
 --wg=WG     Name of Whole Genome file to summarize
 explain     Explain the file format
"""
import gzip
import struct
import logging
logger = logging.getLogger(__name__)

__version__ = '0.1.0'

__explain__ = """

File format

    Header-----------------------------------------------
    [char10]   - version string of wg format
    [char255]  - Human readable species name (string)
    [uint16]   - max number of chromosomes (max 65535)
    [uint16]   - actual number of chromosomes in file

    Index-------------------------------------------------
    For each chromosome copy --------------------------
      [uint16]  - chromosome number (1,2,3,4 ...)
      [uint16]  - chromosome copy (1,2,3 ...)
      [char255] - GenBank accession other standard and unique identifier
      [uint64]  - start byte of data in this file
      [uint64]  - length of sequence

    Data--------------------------------------------------
    For each chromosome copy --------------------------
      [uchar]   - Nucleotide data
        ...

For WholeGenomePos files the Data is uint64 rather than uchar. Everything else is the same.
"""


def cik(chrom_no=1, chrom_cpy=1):
  return '{:d}:{:d}'.format(chrom_no, chrom_cpy)


class WholeGenome():
  """This takes care of reading in an existing WG file or creates an empty one for writing to disk.

  >>> import tempfile
  >>> fname = tempfile.mktemp()
  >>> with WholeGenome(fname, chrom_count=4, compress_level=4) as wg:
  ...   wg.insert_seq('GATTACA', 1,1)
  ...   wg.insert_seq('GATTACA', 1,1)  # Can't do this twice, will be caught
  ...   wg.insert_seq('GATTACAGATTACA', 2,1)
  ...   wg.insert_seq('GATTACT', 1,2)
  ...   wg.insert_seq('GATTACTGA', 2,2)
  ...   wg.insert_seq('GATTACTGA', 3,1)  # Can't do this - too many chromosomes now, will be caught
  True
  False
  True
  True
  True
  False
  >>> with WholeGenome(fname) as wg:
  ...   print wg[1,1]
  ...   print wg['2:2']
  ...   print wg[4,2]  # Will return none - no such chromosome copy
  ('GATTACA', 'Test')
  ('GATTACTGA', 'Test')
  (None, 'No such sequence')
  """
  header_fmt = '10s 255s H H'
  index_fmt = 'H H 255s Q Q'

  def __init__(self, fname=None, species=None, chrom_count=None, compress_level=6):
    """Create a new instance of a whole genome
    Inputs:
      fname          - Open this file for reading if species and chrom_count are None otherwise open for
                       writing and initialize with given chrom_count
      species        - Ascii string corresponding to name of species. Leave None to load
      chrom_count    - Total number of chromosomes in genome. Leave None to load
      """
    if fname is None:  # Create an empty file
      self.header = {
        'version': __version__,
        'species': 'Empty',
        'max chromosome count': 0,
        'actual chromosome count': 0
      }
      self.index = {}
      self.writing = False
      return

    # The only tricky thing here is that when writing we can't use a gzipped file (we can't rewind). So, for writing
    # we open a regular file and we gzip it on __del__ or close()
    mode = 'wb' if chrom_count else 'rb'
    if mode == 'wb':  # Need to write a new file
      self.fp = open(fname + '.tmp', 'w+b')  # open(fname, mode='wb')
      assert compress_level > 0
      assert compress_level < 10
      self.compress_level = int(compress_level)  # We'll compress it on exit
      self.fname = fname
      self.writing = True
    else:  # Opening existing file (gzipped)
      self.fp = gzip.open(fname, mode='rb', compresslevel=compress_level or 9)
      self.writing = False

    if mode == 'rb':  # If this is an existing file we should load the header and indexes
      self.header = self.read_header()
      self.index = self.read_index()
      self.reverse_index = self.create_reverse_index()
    else:  # We should initialize the file
      self.header = {
        'version': __version__,
        'species': species or 'Test',
        'max chromosome count': chrom_count,
        'actual chromosome count': 0
      }
      self.index = {}
      self.write_header()
      # Some pointer locations that are important when we write to a file
      self.index_pos = struct.calcsize(self.header_fmt)  # Our index starts from here
      self.seq_data_pos = struct.calcsize(self.index_fmt) * self.header['max chromosome count'] + self.index_pos

  def close(self):
    if self.writing:
      self.fp.close()  # Flush data
      import subprocess
      import os
      subprocess.call(['gzip', self.fp.name, '-{:d}'.format(self.compress_level)])  # Use the right tool for the job
      os.rename(self.fp.name + '.gz', self.fname)

  # __enter__ and __exit__ are needed for Python contexts ('with')
  def __enter__(self):
    return self

  def __exit__(self, type, value, traceback):
    self.close()
    return True if type is None else False

  def __len__(self):
    return len(self.index)

  def read_header(self):
    try:
      self.fp.seek(0)
    except AttributeError:
      logger.error("Can't read header: File not initialized")
      return None

    values = struct.unpack(self.header_fmt, self.fp.read(struct.calcsize(self.header_fmt)))
    return {
      'version': values[0].strip('\0'),
      'species': values[1].strip('\0'),
      'max chromosome count': values[2],
      'actual chromosome count': values[3]
    }

  def write_header(self):
    self.fp.seek(0)
    self.fp.write(struct.pack(self.header_fmt, self.header['version'], self.header['species'][:255],
                              self.header['max chromosome count'], self.header['actual chromosome count']))

  def read_index(self):
    def read_index_entry(fp):
      values = struct.unpack(self.index_fmt, fp.read(struct.calcsize(self.index_fmt)))
      return {
        'chromosome number': values[0],
        'chromosome copy': values[1],
        'sequence id': values[2].strip('\0'),
        'start byte of sequence data': values[3],
        'length of sequence': values[4]
      }

    try:
      self.fp.seek(struct.calcsize(self.header_fmt))
    except AttributeError:
      logger.error("Can't read index: File not initialized")
      return None

    index = {}
    for n in range(self.header['actual chromosome count']):
      this_index = read_index_entry(self.fp)
      index['{:d}:{:d}'.format(this_index['chromosome number'], this_index['chromosome copy'])] = this_index
    return index

  def create_reverse_index(self):
    """When aligners work on the reads for our data they return data labeled using the first part of the sequence id
    In order to efficiently find which chromosome key a sequence id refers to we create this reverse index when we
    load the data."""
    return {v['sequence id'].split()[0]: tuple([int(m) for m in k.split(':')]) for k, v in self.index.iteritems()}

  def sorted_index_keys(self):
    """Return the index (chromosome) keys sorted by number."""
    return sorted(self.index, key=lambda k: int(k.split(':')[0]))

  def append_index(self, chrom_no=1, chrom_cpy=1, seq_id='Test', start_byte=0, seq_len=0):
    self.fp.seek(self.index_pos)
    self.fp.write(struct.pack(self.index_fmt, chrom_no, chrom_cpy, seq_id, start_byte, seq_len))
    self.index_pos = self.fp.tell()
    self.header['actual chromosome count'] += 1
    self.write_header()
    self.index = self.read_index() # We mirror this for internal book-keeping (keep track of what has been added and how many)

  def insert_seq(self, seq, chrom_no=1, chrom_cpy=1, seq_id='Test'):
    """Given a sequence, insert this into the file and update the index."""
    # Check if this already exists in the index. If so, return an error.
    if cik(chrom_no, chrom_cpy) in self.index:
      logger.error('This sequence has already been inserted')
      return False

    if len(self.index) == self.header['max chromosome count']:
      logger.error('This file is set for {:d} chromosomes and that count has been already reached'.format(self.header['max chromosome count']))
      return False

    self.fp.seek(self.seq_data_pos)
    start_byte = self.fp.tell()
    self._mywrite(seq)
    self.seq_data_pos = self.fp.tell()

    self.append_index(chrom_no, chrom_cpy, seq_id, start_byte, len(seq))

    return True

  def _mywrite(self, seq):
    self.fp.write(seq)

  def __getitem__(self, item):
    """Cute ways to get our sequence data."""
    key = cik(item[0], item[1]) if isinstance(item, tuple) else item
    if key in self.index:
      this_index = self.index[key]
      self.fp.seek(this_index['start byte of sequence data'])
      return self._myread(this_index['length of sequence']), this_index['sequence id']
    else:
      logger.error('No such sequence {:s}'.format(key))
      return None, 'No such sequence'

  def _myread(self, n):
    return self.fp.read(n)


class WholeGenomePos(WholeGenome):
  """This class is used to store POS information for mutated genomes to enable reads.py to produce proper CIGARS.

  >>> import tempfile
  >>> fname = tempfile.mktemp()
  >>> with WholeGenomePos(fname, chrom_count=4) as wg:
  ...   wg.insert_seq([0,1,2,3,4,5,6], 1,1)
  ...   wg.insert_seq([0,1,2,3], 2,1)
  True
  True
  >>> with WholeGenomePos(fname) as wg:
  ...   print wg[1,1]
  ...   print wg['2:1']
  ...   print wg[4,2]  # Will return none - no such chromosome copy
  ((0, 1, 2, 3, 4, 5, 6), 'Test')
  ((0, 1, 2, 3), 'Test')
  (None, 'No such sequence')
  """
  def _mywrite(self, seq):
    self.fp.write(struct.pack('I' * len(seq), *seq))

  def _myread(self, n):
    return struct.unpack('I' * n, self.fp.read(4 * n))


# class WholeGenomeStruct(WholeGenome):
#   """A generalization of WholeGenome to store any kind of data.
#
#   We really should store the 'fmt' string in the header of the file. This will break compatibility with prev version of
#   file, so will do later
#
#   >>> import tempfile
#   >>> fname = tempfile.mktemp()
#   >>> with WholeGenomeStruct(fname, fmt='f I', chrom_count=2) as wg:
#   ...   wg.insert_seq([(0.0, 0), (0.1, 1), (0.2, 2), (0.3, 3), (0.4, 4), (0.5, 5), (0.6, 6)], 1,1)
#   ...   wg.insert_seq([(0.0, 0), (0.6, 6)], 2,1)
#   True
#   True
#   >>> with WholeGenomeStruct(fname, fmt='f I') as wg:
#   ...   print wg[1,1]
#   ...   print wg['2:1']
#   ...   print wg[4,2]  # Will return none - no such chromosome copy
#   ((0, 1, 2, 3, 4, 5, 6), 'Test')
#   ((0, 1, 2, 3), 'Test')
#   (None, 'No such sequence')
#   """
#   import itertools
#
#   def __init__(self, *args, **kwargs):
#     self.fmt = kwargs.pop('fmt')
#     self.fmt_size = struct.calcsize(self.fmt)
#     WholeGenome.__init__(self, *args, **kwargs)
#
#   def _mywrite(self, seq):
#     self.fp.write(struct.pack(self.fmt * len(seq), *itertools.chain.from_iterable(seq)))
#
#   def _myread(self, n):
#     return struct.unpack(self.fmt * n, self.fp.read(self.fmt_size * n))

if __name__ == "__main__":
  import json
  import docopt
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  if args['explain']:
    print __explain__
    exit(0)

  with WholeGenome(fname=args['--wg']) as wg:
    print 'Header'
    print json.dumps(wg.header, indent=2)
    print 'Index'
    print json.dumps(wg.index, indent=2)
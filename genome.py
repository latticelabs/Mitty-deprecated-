"""This module defines a whole genome file (.wg) that allows us to store multiple sequences in an organized fashion. It
also defines an interface to read in chosen sequences from the file. The file is meant to represent an individual and
carries multiple chromosomes and multiple copies of each chromosome. """
import tempfile
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
    [uint16]   - number of chromosomes max 65535

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
"""


def cik(chrom_no=1, chrom_cpy=1):
  return '{:d}:{:d}'.format(chrom_no, chrom_cpy)


class WholeGenome():
  """This takes care of reading in an existing WG file or creates an empty one for writing to disk.

  >>> import tempfile
  >>> fname = tempfile.mktemp()
  >>> with WholeGenome(fname, chrom_count=4) as wg:
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
  ...   print wg.get_seq(1,1)
  ...   print wg.get_seq(2,2)
  ...   print wg.get_seq(4,2)  # Will return none - no such chromosome copy
  GATTACA
  GATTACTGA
  None
  """
  header_fmt = '10s 255s H'
  index_fmt = 'H H 255s Q Q'

  def __init__(self, fname=None, species=None, chrom_count=None, compress_level=None):
    """Create a new instance of a whole genome
    Inputs:
      fname          - Open this file for reading if species and chrom_count are None otherwise open for
                       writing and initialize with given chrom_count
      species        - Ascii string corresponding to name of species. Leave None to load
      chrom_count    - Total number of chromosomes in genome. Leave None to load
      """
    if fname is None:
      return
    # The only tricky thing here is that when writing we can't use a gzipped file (we can't rewind). So, for writing
    # we open a regular file and we gzip it on __del__ or close()
    mode = 'wb' if chrom_count else 'rb'
    if mode == 'wb':  # Need to write a new file
      self.fp = tempfile.TemporaryFile()  #  open(fname, mode='wb')
      self.compress_level = compress_level  # We'll compress it on exit
      self.fname = fname
      self.writing = True
    else:  # Opening existing file (gzipped)
      self.fp = gzip.open(fname, mode='rb', compresslevel=compress_level or 9)
      self.writing = False

    if mode == 'rb':  # If this is an existing file we should load the header and indexes
      self.header = self.read_header()
      self.index = self.read_index()
    else:  # We should initialize the file
      self.header = {
        'version': __version__,
        'species': species or 'Test',
        'chromosome count': chrom_count
      }
      self.index = {}
      self.write_header()
      # Some pointer locations that are important when we write to a file
      self.index_pos = struct.calcsize(self.header_fmt)  # Our index starts from here
      self.seq_data_pos = struct.calcsize(self.index_fmt) * self.header['chromosome count'] + self.index_pos

  def __del__(self):
    """Only needed if we are writing. We'd like to compress on exit."""
    self.close()

  def close(self):
    if self.writing:
      block_size = 1073741824
      with gzip.open(self.fname, mode='wb', compresslevel=self.compress_level or 9) as fp:
        self.fp.seek(0)
        bytes_written = 1
        while bytes_written:
          bytes_written = fp.write(self.fp.read(block_size))  # This should compress the file.

  # __enter__ and __exit__ are needed for Python contexts ('with')
  def __enter__(self):
    return self

  def __exit__(self, type, value, traceback):
    self.close()

  def read_header(self):
    try:
      self.fp.seek(0)
    except AttributeError:
      logger.error("Can't read header: File not initialized")
      return None

    values = struct.unpack(self.header_fmt, self.fp.read(struct.calcsize(self.header_fmt)))
    return {
      'version': values[0],
      'species': values[1],
      'chromosome count': values[2]
    }

  def write_header(self):
    self.fp.seek(0)
    self.fp.write(struct.pack(self.header_fmt, self.header['version'],
                              self.header['species'], self.header['chromosome count']))

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
    for n in range(self.header['chromosome count']):
      this_index = read_index_entry(self.fp)
      index['{:d}:{:d}'.format(this_index['chromosome number'], this_index['chromosome copy'])] = this_index
    return index

  def append_index(self, chrom_no=1, chrom_cpy=1, seq_id='Test', start_byte=0, seq_len=0):
    self.fp.seek(self.index_pos)
    self.fp.write(struct.pack(self.index_fmt, chrom_no, chrom_cpy, seq_id, start_byte, seq_len))
    self.index_pos = self.fp.tell()
    self.index[cik(chrom_no, chrom_cpy)] = {
        'chromosome number': chrom_no,
        'chromosome copy': chrom_cpy,
        'sequence id': seq_id,
        'start byte of sequence data': start_byte,
        'length of sequence': seq_len
    }  # We mirror this for internal book-keeping (keep track of what has been added and how many)

  def insert_seq(self, seq, chrom_no=1, chrom_cpy=1, seq_id='Test'):
    """Given a sequence, insert this into the file and update the index."""
    # Check if this already exists in the index. If so, return an error.
    if cik(chrom_no, chrom_cpy) in self.index:
      logger.error('This sequence has already been inserted')
      return False

    if len(self.index) == self.header['chromosome count']:
      logger.error('This file is set for {:d} chromosomes and that count has been already reached'.format(self.header['chromosome count']))
      return False

    self.fp.seek(self.seq_data_pos)
    start_byte = self.fp.tell()
    self.fp.write(seq)
    self.seq_data_pos = self.fp.tell()

    self.append_index(chrom_no, chrom_cpy, seq_id, start_byte, len(seq))

    return True

  def get_seq(self, chrom_no=1, chrom_cpy=1):
    if cik(chrom_no, chrom_cpy) in self.index:
      self.fp.seek(self.index[cik(chrom_no, chrom_cpy)]['start byte of sequence data'])
      return self.fp.read(self.index[cik(chrom_no, chrom_cpy)]['length of sequence'])
    else:
      logger.error('No such sequence {:s}'.format(cik(chrom_no, chrom_cpy)))
      return None
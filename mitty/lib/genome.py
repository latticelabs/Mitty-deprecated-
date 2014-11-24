"""This module contains a class that represents a reference genome and provides us with an abstraction to handle input
genome data"""
import numpy
import glob
import re
import logging
logger = logging.getLogger(__name__)


def split_multi_fasta_gz(fa_fname, dir_out):
  """Given a gzipped multi fa.gz file split it into separate files as used by mitty"""
  def write_it_out(d_out, ch, sid, seq):
    with open(os.path.join(d_out, 'chr{:d}.fa'.format(ch)), 'w') as fp_out:
      fp_out.write(sid + '\n')
      fp_out.writelines(seq)

  import gzip
  import os
  with gzip.open(fa_fname, 'rb') as fp:
    this_seq = []
    chrom = 0
    for line in fp:
      line = line.strip()
      if len(line) == 0: continue
      if line[0] == '>':
        if len(this_seq):
          write_it_out(dir_out, chrom, seq_id, this_seq)
        seq_id = line[1:]
        this_seq = []
        chrom += 1
      else:
        this_seq += [line.upper()]
    if len(this_seq):
      write_it_out(dir_out, chrom, seq_id, this_seq)


def convert_fasta(fa_fname):
  """Read a fasta with newlines and lowercase letters. Strip out newlines and convert all to uppercase."""
  with open(fa_fname, 'r') as fasta_fp:
    seq_id = fasta_fp.readline()
    lines = fasta_fp.readlines()

  with open(fa_fname, 'w') as fasta_fp:
    fasta_fp.write(seq_id)
    for line in lines:
      fasta_fp.write(line.strip().upper())


# Discovered that reading an unzipped file is WAAY faster that reading a zipped one
def load_single_fasta(fa_fname, as_numpy=False):
  """Expects a fasta file with only one sequence and only upper case letters - will read other files but the result
  is not sanitized in any way - newlines and repeat masks are left in.
  if as_numpy is set we will get the result as a numpy char array"""
  with open(fa_fname, 'r') as fasta_fp:
    seq_id = fasta_fp.readline()[1:-1]
    logger.debug('Read {:s}'.format(seq_id))
    seq = numpy.fromfile(fasta_fp, dtype='u1') if as_numpy else fasta_fp.read()
  return seq


class FastaGenome():
  """This represents a genome as a collection of fasta files stored in a directory.

  The fasta files should be numbered chr1, chr2, chr3 ...
  Sequences can be accessed using the slice notation e.g.::

    ref = FastaGenome('examples/data', persist=True)
    chr1 = ref[1]  # Loads the sequence data and keeps it attached to the object

  """
  def __init__(self, seq_dir='', persist=False):
    """
    :param str seq_dir: path to the directory containing the fasta files
    :param bool persist: If true, once loaded, a sequence is kept in memory. This speeds up operations at the
                           expense of memory."""
    self.dir = seq_dir
    self.persist = persist
    self.sequences = {}

  def __getitem__(self, item):
    def load_seq(it):
      fa_fname = glob.os.path.join(self.dir, 'chr{:s}.fa'.format(str(item)))
      if glob.os.path.exists(fa_fname):
        sq = load_single_fasta(fa_fname)
      else:
        logger.warning('{:s} does not exist'.format(fa_fname))
        sq = None
      return sq

    if self.persist:
      if item not in self.sequences:
        seq = load_seq(item)
        self.sequences[item] = seq
      else:
        seq = self.sequences[item]
    else:
      seq = load_seq(item)
    return seq

  def __len__(self):
    return len(self.genome_header())

  def get(self, item, default=None):
    return self.__getitem__(item) or default

  def sorted_chrom_idx(self):
    """Get a list of integer indexes of files in the directory of the format chrX where X is an integer."""
    def f_idx(fn0):
      m = re.search(r'chr([0-9]+).fa$', fn0)
      return int(m.groups()[0]) if m else None
    return sorted(filter(lambda x: x is not None, [f_idx(fn) for fn in glob.os.listdir(self.dir)]))

  def genome_header(self):
    """Return a list of tuples of the form (seq_id, seq_len, seq_offset) corresponding to the files in the directory.
    Files not in the format chrX where X is an integer are ignored"""
    def process_file(fname, seq_offset=[0]):
      # http://stackoverflow.com/questions/3432830/list-comprehension-for-running-total
      try:
        with open(fname, 'r') as fasta_fp:
          seq_id = fasta_fp.readline()[1:-1]
          idx = fasta_fp.tell()
          fasta_fp.seek(0, 2)
          seq_len = fasta_fp.tell() - idx
          this_seq_offset = seq_offset[0]
          seq_offset[0] += seq_len
      except IOError:
        seq_id, seq_len, this_seq_offset = None, None, None
      return seq_id, seq_len, this_seq_offset

    return filter(lambda x: x[0] is not None, [process_file(glob.os.path.join(self.dir, 'chr{:s}.fa'.format(str(chrom)))) for chrom in self.sorted_chrom_idx()])
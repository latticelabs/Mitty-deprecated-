"""This module contains a class that represents a reference genome and provides us with an abstraction to handle input
genome data"""
import numpy
import glob
import logging
logger = logging.getLogger(__name__)


human_chromosomes = list(range(1, 25))


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
  """This represents a genome as a collection of fasta files stored in a directory. The fasta files should be numbered
  chr1, chr2, chr3 ... chr23, chr24. Sequences can be accessed using the slice notation e.g.

  ref = FastaGenome('examples/data')
  ch = ref[1]  # Loads the data

  """
  def __init__(self, seq_dir='', persist=False):
    """If persist is set we retain a loaded sequence in memory."""
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

  def get(self, item, default=None):
    return self.__getitem__(item) or default

  def genome_header(self):
    """Return a dictionary keyed to all the available chromosomes."""
    def process_sequence(fname):
      try:
        with open(fname, 'r') as fasta_fp:
          seq_id = fasta_fp.readline()[1:-1]
          idx = fasta_fp.tell()
          fasta_fp.seek(0, 2)
          seq_len = fasta_fp.tell() - idx
      except IOError:
        seq_id, seq_len = None, None
      return seq_id, seq_len

    return filter(lambda x: x[0] is not None, [process_sequence(glob.os.path.join(self.dir, 'chr{:s}.fa'.format(str(chrom)))) for chrom in human_chromosomes])
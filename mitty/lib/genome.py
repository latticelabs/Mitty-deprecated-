"""This module contains a class that represents a reference genome and provides us with an abstraction to handle input
genome data"""
import gzip
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
def load_single_fasta(fa_fname):
  """Expects a fasta file with only one sequence and only upper case letters - will read other files but the result
  is not sanitized in any way - newlines and repeat masks are left in."""
  with open(fa_fname, 'r') as fasta_fp:
    seq_id = fasta_fp.readline()[1:-1]
    seq = fasta_fp.read()
  return seq


class FastaGenome():
  """This represents a genome as a collection of fasta files stored in a directory. The fasta files should be numbered
  chr1, chr2, chr3 ... chrX, chrY"""
  def __init__(self, seq_dir='', load_at_once=False, chrom_list=None):
    self.dir = seq_dir
    self.load_at_once = load_at_once
    if self.load_at_once:
      self.sequences = {
        str(chrom): load_single_fasta(glob.os.path.join(self.dir,'chr{:s}.fa'.format(str(chrom))))
        for chrom in chrom_list
      }

  def __getitem__(self, item):
    if self.load_at_once:
      return self.sequences.get(str(item), None)
    else:
      fa_fname = glob.os.path.join(self.dir, 'chr{:s}.fa'.format(str(item)))
      if glob.os.path.exists(fa_fname):
        return load_single_fasta(fa_fname)
      else:
        logger.warning('{:s} does not exist'.format(fa_fname))
        return None

  def get(self, item, default):
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




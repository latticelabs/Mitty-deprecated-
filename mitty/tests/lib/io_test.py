import tempfile
import gzip

import mitty.lib.io as mio
import mitty.lib.variants as vr
from mitty.tests import *  # To get definitions from the setup script
from nose.tools import assert_raises


def unzipped_multi_fasta_test():
  """Load unzipped multi-fasta."""
  ref = mio.Fasta(multi_fasta=os.path.join(example_data_dir, 'chimera.fa'))
  assert len(ref) == 0
  assert len(ref[4]) == 702
  assert len(ref) == 4


def gzipped_multi_fasta_test():
  """Load gzipped multi-fasta."""
  ref = mio.Fasta(multi_fasta=os.path.join(example_data_dir, 'chimera.fa.gz'))
  assert len(ref) == 0
  assert len(ref[4]) == 702
  assert len(ref) == 4


def multi_dir_test():
  """Load reference from directory"""
  ref = mio.Fasta(multi_dir=example_data_dir)
  assert len(ref) == 0
  assert len(ref[4]) == 702
  assert len(ref) == 1
  assert len(ref[4]) == 702
  assert len(ref) == 1
  assert len(ref[3]) == 717
  assert len(ref) == 2


class Variant:
  """A lightweight test class that represents a row of a VCF file."""
  def __init__(self, chrom, pos, vid, ref, alt, qual, filter, huh, huh2, gt):
    self.chrom, self.pos, self.ref, self.alt, self.gt = chrom, int(pos), ref, alt, gt.strip()


def simple_vcf_reader(fname):
  """The current version only needs pyvcf for this one test, so we got rid of it, at least for this version"""
  with gzip.open(fname, 'r') as fp:
    # We simply skip the headers
    for line in fp:
      if not line.startswith('#'):
        break
    variants = [Variant(*line.split('\t'))]
    for line in fp:
      variants += [Variant(*line.split('\t'))]
  return variants


def vcf_contextmanager_test():
  """VCF context manager/writing"""
  def do_this(tn):
    with mio.vcf_for_writing(tn, ['a', 'b']) as fp:
      fp.write('a')

  temp_fp, temp_name = tempfile.mkstemp(suffix='.vcf.gz')
  os.close(temp_fp)

  assert_raises(NotImplementedError, do_this, temp_name)

  pos = [1, 10, 20, 30]
  stop = [2, 11, 21, 35]
  ref = ['A', 'C', 'T', 'GAAAA']
  alt = ['AA', 'CAT', 'G', 'G']
  p = [0.1, 0.5, 0.9, 0.2]
  ml = vr.VariantList(pos, stop, ref, alt, p)
  chrom = [(0, 0), (2, 1), (3, 2)]

  with mio.vcf_for_writing(temp_name, ['a']) as fp:
    mio.write_chromosomes_to_vcf(fp, seq_id='chr2', chrom_list=[chrom], master_list=ml)

  v = simple_vcf_reader(temp_name)

  assert v[0].pos == 2
  assert v[0].ref == 'A'
  assert v[0].alt == 'AA'
  assert v[0].gt == '1|0'

  assert v[2].pos == 31
  assert v[2].ref == 'GAAAA'
  assert v[2].alt == 'G'
  assert v[2].gt == '1|1'
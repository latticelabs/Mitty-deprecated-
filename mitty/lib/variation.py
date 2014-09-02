"""This module defines a single, simple, named_tuple that is used in several places by Mitty and should be used by
all the mutation plugins

(start, stop, REF, ALT, het)

"""
from os.path import splitext
from ctypes import *
import pysam
import logging
logger = logging.getLogger(__name__)

# Types of het
HOMOZYGOUS = 0
HET1 = 1
HET2 = 2
GT = ['1/1', '1/0', '0/1']  # This needs to match rev 3 definitions


class Variation(Structure):
  _fields_ = [("POS", c_int32),
              ("stop", c_int32),
              ("REF", c_char_p),
              ("ALT", c_char_p),
              ("het", c_uint8)]

  def __eq__(self, other):
    # This does not do a isinstance check for speed reasons.
    if self.POS == other.POS and self.stop == other.stop and self.REF == other.REF and self.ALT == other.ALT and self.het == other.het:
      return True
    else:
      return False

  def __ne__(self, other):
      return not self.__eq__(other)

  def __repr__(self):
    return '(POS={0},stop={1},REF={2},ALT={3},het={4})'.format(self.POS, self.stop, self.REF, self.ALT, GT[self.het])


def vcopy(v, het=None):
  if het is None:
    return Variation(v.POS, v.stop, v.REF, v.ALT, v.het)
  else:
    return Variation(v.POS, v.stop, v.REF, v.ALT, het)


def compress_and_index_vcf(in_vcf_name, out_vcf_name):
  #bgzip -c sorted.vcf > sorted.vcf.gz
  #tabix sorted.vcf.gz
  logger.debug('Compressing and indexing {:s} to {:s}'.format(in_vcf_name, out_vcf_name))
  pysam.tabix_compress(in_vcf_name, out_vcf_name, force=True)
  pysam.tabix_index(out_vcf_name, force=True, preset='vcf')


def vcf2chrom(vcf_rdr):
  """Given a vcf reader corresponding to one chromosome, read in the variant descriptions into our format. The result is
  sorted if the vcf file is sorted.
  """
  chrom = []
  append = chrom.append
  for variant in vcf_rdr:
    alt = variant.ALT[0].sequence if variant.ALT[0] is not None else ''
    ref = variant.REF or ''
    start = variant.POS  # Note, we are in VCF coordinates!
    stop = variant.POS + len(ref)
    het = HOMOZYGOUS

    try:
      if variant.samples[0].gt_nums[0] == '0':
        het = HET2
      if variant.samples[0].gt_nums[2] == '0':
        if het == HET2:  # 0/0 means this does not exist in this sample
          continue
        else:
          het = HET1
    except IndexError:  # No genotype info, will assume homozygous
        pass

    append(Variation(start, stop, ref, alt, het))

  return chrom


def parse_vcf(vcf_rdr, chrom_list):
  """Given a vcf reader load in all the chromosomes."""
  g1 = {}
  for chrom in chrom_list:
    try:
      g1[chrom] = vcf2chrom(vcf_rdr.fetch(chrom, start=0))
    except KeyError:
      g1[chrom] = []

  return g1


def vcf_save(g1, fp, sample_name='sample'):
  """Given a genome save it to a VCF file."""
  # Write header
  fp.write(
    "##fileformat=VCFv4.1\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{:s}\n".format(sample_name)
  )
  # Write lines
  wr = fp.write
  for chrom, variants in g1.iteritems():  # We don't bother sorting - we do that when we bgzip and index it
    ch = str(chrom)
    for var in variants:
      # In the VCF file no REF or ALT is indicated by a .
      ref = var.REF if var.REF != '' else '.'
      alt = var.ALT if var.ALT != '' else '.'
      #  CHROM    POS   ID   REF   ALT   QUAL FILTER INFO FORMAT tsample
      wr(ch + "\t" + str(var.POS) + "\t.\t" + ref + "\t" + alt + "\t100\tPASS\t.\tGT\t" + GT[var.het] + "\n")


def vcf_save_gz(g1, vcf_gz_name, sample_name='sample'):
  """Save .vcf, bgzip and index it. File name should have .gz at the end, but it's not a drama if doesnt. Sigh"""
  vcf_name, ext = splitext(vcf_gz_name)
  if ext != '.gz':  # Like I said, not a drama
    vcf_name += '_srt.vcf'
    vcf_gz_name += '.gz'

  with open(vcf_name, 'w') as fp:
    vcf_save(g1, fp, sample_name=sample_name)

  compress_and_index_vcf(vcf_name, vcf_gz_name)
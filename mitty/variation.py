"""This module defines a single, simple, named_tuple that is used in several places by Mitty and should be used by
all the mutation plugins

(start, stop, REF, ALT, het)

"""
import vcf
import pysam
import subprocess
from collections import namedtuple
import logging
logger = logging.getLogger(__name__)
# Types of het
HOMOZYGOUS = 0
HET1 = 1
HET2 = 2
GT = ['1/1', '1/0', '0/1']  # This needs to match rev 3 definitions
Variation = namedtuple('Variation', 'POS, stop, REF, ALT, het')


def sort_vcf(in_vcf_name, out_vcf_name):
  #vcf-sort the.vcf > sorted.vcf
  logger.debug('Sorting {:s}'.format(in_vcf_name))
  with open(out_vcf_name, 'w') as fp:
    subprocess.call(['vcf-sort', in_vcf_name], stdout=fp)


def compress_and_index_vcf(in_vcf_name, out_vcf_name):
  #bgzip -c sorted.vcf > sorted.vcf.gz
  #tabix sorted.vcf.gz
  logger.debug('Compressing and indexing {:s}'.format(in_vcf_name))
  pysam.tabix_compress(in_vcf_name, out_vcf_name, force=True)
  pysam.tabix_index(out_vcf_name, force=True, preset='vcf')


def vcf2chrom(vcf_rdr):
  """Given a vcf reader corresponding to one chromosome, read in the variant descriptions into our format. The result is
  sorted if the vcf file is sorted.
  """
  chrom = []

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

    chrom += [Variation(start, stop, ref, alt, het)]

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


def vcf_save(g1, fp):
  """Given a genome save it to a VCF file."""
  # Write header
  fp.write(
    "##fileformat=VCFv4.1\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
  )
  # Write lines
  for chrom, variants in g1.iteritems():  # We don't bother sorting - we do that when we bgzip and index it
    for var in variants:
      # In the VCF file no REF or ALT is indicated by a .
      var_ = var._replace(REF=var.REF or '.', ALT=var.ALT or '.')
      #  CHROM    POS   ID   REF   ALT   QUAL FILTER INFO FORMAT tsample
      fp.write("{chrom}\t{POS}\t.\t{REF}\t{ALT}\t100\tPASS\t.\tGT\t{gt}\n".format(chrom=chrom, gt=GT[var.het], **var_.__dict__))


def vcf_save_gz(g1, vcf_gz_name):
  """Also sort it, bgzip and index it. File name should have .gz at the end, but it's not a drama if doesnt. Sigh"""
  import tempfile, os

  vcf_name, ext = os.path.splitext(vcf_gz_name)
  if ext != '.gz':  # Like I said, not a drama
    vcf_name += '_srt.vcf'

  temp_vcf_fp, temp_vcf_name = tempfile.mkstemp(suffix='.vcf')
  vcf_save(g1, temp_vcf_fp)

  #with open(temp_vcf_name, 'w') as fp:
  #  vcf_save(g1, fp)

  sort_vcf(temp_vcf_name, vcf_name)
  compress_and_index_vcf(vcf_name, vcf_gz_name)
  os.remove(temp_vcf_name)
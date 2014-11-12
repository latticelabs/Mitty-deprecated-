"""This module defines a single, simple, named_tuple that is used in several places by Mitty and should be used by
all the mutation plugins

(start, stop, REF, ALT, het)

"""
from os.path import splitext
#from ctypes import *
import pysam
import logging
logger = logging.getLogger(__name__)

# Types of het
HOMOZYGOUS = 0
HET1 = 1
HET2 = 2
GT = ['1/1', '1/0', '0/1']  # This needs to match rev 3 definitions



cdef class Variation:
  cdef public:
    int POS, stop
    str REF, ALT
    char het

  def __cinit__(self, int _pos=0, int _stop=0, str _ref='', str _alt='', char _het=0):
    self.POS, self.stop, self.REF, self.ALT, self.het = _pos, _stop, _ref, _alt, _het

  def __richcmp__(self, Variation other, int op):
    if op == 2:
      return self.POS == other.POS and self.stop == other.stop and self.REF == other.REF and self.ALT == other.ALT and self.het == other.het
    elif op == 3:
      return not self.__richcmp__(other, 2)
    else:
      return False

  def __repr__(self):
    return '(POS={0},stop={1},REF={2},ALT={3},het={4})'.format(self.POS, self.stop, self.REF, self.ALT, GT[self.het])


def copy_genome(g1):
  """g1 - dictionary with chromosome name as key, each value is list of Variations."""
  return {k: copy_chromosome(v) for k, v in g1.iteritems()}


def copy_chromosome(c1):
  """c1 - list of Variations corresponding to chromosome."""
  return [vcopy(v) for v in c1]


cpdef inline Variation vcopy(Variation v, het=None):
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
    except (ValueError, KeyError):  # New version of pyvcf changed the error
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
    vcf_name += '.vcf'
    vcf_gz_name = vcf_name + '.gz'

  with open(vcf_name, 'w') as fp:
    vcf_save(g1, fp, sample_name=sample_name)

  compress_and_index_vcf(str(vcf_name), str(vcf_gz_name))
  # tabix can't understand unicode, needs bytes


cdef inline bint overlap(Variation x, Variation y):
  if x is None or y is None: return False
  if y.POS - 1 <= x.POS <= y.stop + 1 or y.POS - 1 <= x.stop <= y.stop + 1 or x.POS <= y.POS - 1 <= y.stop + 1 <= x.stop:
    # Potential overlap
    if x.het == y.het or x.het == HOMOZYGOUS or y.het == HOMOZYGOUS:
      return True
  return False


def merge_variants_with_chromosome(c1, dnv):
  """
  Given an exiting chromosome (in variant format) merge new variants into it in zipper fashion

  Args:
    c1 (variant list): The original chromosome
    dnv (variant list): The proposed variants

  Returns:
    c2 (variant list): The resultant chromosome with variant collisions arbitrated

  Algorithm:
    o(x,y) = True if x and y overlap
           = False otherwise
    e = existing list of variants
    d = denovo list of variants
    n = new list of variants being built

    o(e, d) = True: add(e) e++, d++
            = False:
              e < d ? add(e) e++
              else:
                o(n, d) = True: d++
                        = False: add(d), d++

  """
  c1_iter, dnv_iter = c1.__iter__(), dnv.__iter__()
  c2 = []
  append = c2.append
  last_new = None
  # Try the zipper
  cdef Variation existing = next(c1_iter, None), denovo = next(dnv_iter, None)
  while existing is not None and denovo is not None:
    if overlap(existing, denovo):
      # This will collide, resolve in favor of existing and advance both lists
      append(vcopy(existing))
      last_new = existing
      existing, denovo = next(c1_iter, None), next(dnv_iter, None)
    else:
      if existing.POS <= denovo.POS:  # Zip-in existing
        append(vcopy(existing))
        last_new = existing
        existing = next(c1_iter, None)
      else:  # Can we zip-in denovo?
        if not overlap(last_new, denovo):
          append(vcopy(denovo))
          last_new = denovo
        denovo = next(dnv_iter, None)  # In either case, we need to advance denovo

  # Now pick up any slack
  if existing is not None:  # Smooth sailing, just copy over the rest
    while existing is not None:
      append(vcopy(existing))
      existing = next(c1_iter, None)
  else:  # Need to test for overlap before copying over
    while denovo is not None:
      if not overlap(last_new, denovo):
        append(vcopy(denovo))
        last_new = c2[-1]
      denovo = next(dnv_iter, None)  # In either case, we need to advance denovo

  return c2
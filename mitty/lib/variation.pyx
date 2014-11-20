"""This module defines a structure to carry variation information and some functions to perform operations on
genomes (collections of variations). For a description of design choices and algorithms please see the developer
documentation
"""
from os.path import splitext
import pysam
import logging
logger = logging.getLogger(__name__)

# Types of het
ABSENT = 0
HOMOZYGOUS = 3
HET_01 = 1
HET_10 = 2
#       00     01     10     11
GT = ['0/0', '0/1', '1/0', '1/1']  # This needs to match rev 3 definitions


cdef class VariationData:
  """Carries POS, stop REF and ALT information."""
  cdef public:
    long int POS, stop
    str REF, ALT

  def __cinit__(self, int _pos=0, int _stop=0, str _ref='', str _alt=''):
    self.POS, self.stop, self.REF, self.ALT = _pos, _stop, _ref, _alt


cdef inline char encode_variation_metadata(char het=0, char recessive=0, float fitness=0):
  return (het << 6) | (recessive << 5) | (<char>(fitness * 15 + 16) & 0x1f)


cdef inline char decode_het_info(char metadata):
  return metadata >> 6


def decode_variation_metadata(char metadata):
  return metadata >> 6, (metadata >> 5) & 0x1, ((metadata & 0x1f) - 16)/15.0


def new_variation(int _pos=0, int _stop=0, str _ref='', str _alt='', char het=0, char recessive=0, float fitness=0):
  return Variation(VariationData(_pos, _stop, _ref, _alt), encode_variation_metadata(het, recessive, fitness))


cdef class Variation:
  """Carries *reference* to a VariationData instance and heterozygosity, recessiveness and fitness data."""
  cdef public:
    VariationData vd
    char metadata
    # This 8 bit field carries the following information
    # 7,6   - het
    # 5     - recessive
    # 0-4   - fitness value 0-15 -ve 17-31 +ve

  # http://docs.cython.org/src/userguide/special_methods.html
  def __cinit__(self, VariationData _vd, char _v_info):
    self.vd, self.metadata = _vd, _v_info

  # http://docs.cython.org/src/userguide/special_methods.html
  def __richcmp__(self, Variation other, int op):
    if op == 2:
      return self.vd.POS == other.vd.POS and self.vd.stop == other.vd.stop and \
             self.vd.REF == other.vd.REF and self.vd.ALT == other.vd.ALT and \
             decode_het_info(self.metadata) == decode_het_info(other.metadata)
    elif op == 3:
      return not self.__richcmp__(other, 2)
    else:
      return False

  def __repr__(self):
    """Return a human friendly printed representation of the variant."""
    return '(POS={0},stop={1},REF={2},ALT={3},het={4},r={5},fit={6})'.\
      format(self.vd.POS, self.vd.stop, self.vd.REF, self.vd.ALT, GT[self.metadata >> 6], (self.metadata >> 5) & 0x1,
             ((self.metadata & 0x1f) - 16)/15.0)


cdef inline Variation copy_variant(v1):
  return Variation(v1.vd, v1.metadata)


def copy_genome(g1):
  """g1 - dictionary with chromosome name as key, each value is list of Variations."""
  return {k: copy_chromosome(v) for k, v in g1.iteritems()}


def copy_chromosome(c1):
  """c1 - list of Variations corresponding to a chromosome."""
  return [Variation(v, v.metadata) for v in c1]


def compress_and_index_vcf(in_vcf_name, out_vcf_name):
  """Given an uncompressed, but sorted, vcf, compress and index it."""
  #bgzip -c sorted.vcf > sorted.vcf.gz
  #tabix sorted.vcf.gz
  logger.debug('Compressing and indexing {:s} to {:s}'.format(in_vcf_name, out_vcf_name))
  pysam.tabix_compress(in_vcf_name, out_vcf_name, force=True)
  pysam.tabix_index(out_vcf_name, force=True, preset='vcf')


#TODO implement saving of fitness and recessive fields and loading (if present)
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
        het = HET_01
      if variant.samples[0].gt_nums[2] == '0':
        if het == HET_01:  # 0/0 means this does not exist in this sample
          het = ABSENT
          continue
        else:
          het = HET_10
    except IndexError:  # No genotype info, will assume homozygous
        pass

    append(Variation(VariationData(start, stop, ref, alt), encode_variation_metadata(het)))

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
      ref = var.vd.REF if var.vd.REF != '' else '.'
      alt = var.vd.ALT if var.vd.ALT != '' else '.'
      #  CHROM    POS   ID   REF   ALT   QUAL FILTER INFO FORMAT tsample
      wr(ch + "\t" + str(var.vd.POS) + "\t.\t" + ref + "\t" + alt + "\t100\tPASS\t.\tGT\t" + GT[decode_het_info(var.metadata)] + "\n")


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


cdef inline vcopy(Variation x):
  return Variation(x.vd, x.metadata)


cdef inline bint overlap(Variation x, Variation y):
  # Returns true if the footprints of the variations overlap.
  if x is None or y is None: return False
  if y.vd.POS - 1 <= x.vd.POS <= y.vd.stop + 1 or y.vd.POS - 1 <= x.vd.stop <= y.vd.stop + 1 or \
                  x.vd.POS <= y.vd.POS - 1 <= y.vd.stop + 1 <= x.vd.stop:
    # Potential overlap
    x_het, y_het = decode_het_info(x.metadata), decode_het_info(y.metadata)
    if x_het == y_het or x_het == HOMOZYGOUS or y_het == HOMOZYGOUS:
      return True
  return False


def merge_variants(c1, c2):
  """
  Given an existing chromosome (list of variants) merge a new list of variants into it in zipper fashion

  Args:
    c1 (variant list): The original variants
    c2 (variant list): The proposed new variants

  Returns:
    c3 (variant list): The resultant variant list with collisions arbitrated

  Algorithm:
    o(x,y) = True if x and y overlap
           = False otherwise
    c1 = existing list of variants
    c2 = denovo list of variants
    c3 = new list of variants being built

    o(c1, c2) = True: add(c1) c1++, c2++
              = False:
                c1 < c2 ? add(c1) c1++
                else:
                  o(c3, c2) = True: c2++
                            = False: add(c2), c2++
  """
  c1_iter, dnv_iter = c1.__iter__(), c2.__iter__()
  c3 = []
  append = c3.append
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
        last_new = c3[-1]
      denovo = next(dnv_iter, None)  # In either case, we need to advance denovo

  return c3
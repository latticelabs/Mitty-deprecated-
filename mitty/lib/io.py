"""Some utilities related to loading/saving different file formats, compressing and indexing."""
from os.path import splitext
import os
import gzip
import time

import pysam

import logging
logger = logging.getLogger(__name__)


def load_multi_fasta_gz(fa_fname):
  """Given a gzipped multi fa.gz file load it into a dictionary
  :param fa_fname: fasta.gz file with one or more fasta sequences
  :returns ref_seq: a list of tuples of the form (seq_id, seq_len, seq) in the order they are found in the file, which is the
                    same as chromosome number (with a 1 offset, because of 0-indexing)

  Pure Python 2min 9s to load hg38
  Cythonized 2min 6s - since we are mostly in Python native functions, we are at speed limit
  """
  t0 = time.time()
  ref_seq = []
  with gzip.open(fa_fname, 'r') as fp:
    l = fp.read()
  seq_strings = l.split('>')
  for seq_string in seq_strings:
    if seq_string == '':
      continue
    idx = seq_string.find('\n')
    if idx == -1:
      raise RuntimeError('Something wrong with the fasta file {:s}'.format(fa_fname))
    if idx == 0:
      continue  # Empty line, ignore
    seq_id = seq_string[:idx]
    seq = seq_string[idx:].replace('\n', '').upper()
    ref_seq += [(seq_id, len(seq), seq)]
  t1 = time.time()
  logger.debug('Loaded reference in {:f} s'.format(t1 - t0))
  return ref_seq


def load_vcf(rdr, sample_labels=None, max_rows=-1):
  """UNTESTED. SHOULD WORK

  Given a VCF reader from a particular chromosome, create master list and sample lists.
  :param rdr: VCF reader object (see below)
  :param sample_labels: list of sample names. If None, load all the samples in the file
  :returns l: master list of variants (dict)
  :returns samples: list of list of Genotypes

  fname = '/Users/kghose/Data/vcf1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz'
  vf = vcf.Reader(filename=fname)
  rdr = vf.fetch(chrom=22, start=0)
  """
  if sample_labels is None:
    sample_labels = rdr.samples

  l = {}
  samples = [[] for _ in sample_labels]
  for i, v in enumerate(rdr):
    if i == max_rows:  # Works for -1
      break
    for n, sn in enumerate(sample_labels):
      call = v.genotype(sn)
      if call.is_variant:
        g1 = int(call.gt_alleles[0])
        if g1 == 0:
          g1 = int(call.gt_alleles[1])
      else:
        continue
      g1 -= 1
      try:
        alt = v.ALT[g1].sequence
      except AttributeError:
        continue
      het = HOMOZYGOUS
      if call.gt_alleles[0] == '0':
        het = HET_01
      elif call.gt_alleles[1] == '0':
        het = HET_10
      samples[n] += [cgt(l, v.POS, stop=v.POS + len(v.REF), ref=v.REF, alt=alt, het=het)]
  return l, samples


def save_vcf(vcf_name, chrom, ml, samples, sample_names, bgzip_and_index=True):
  """NOT IMPLEMENTED YET

  Given a master list for a given chromosome and a set of samples, save them to a VCF file.
  :param vcf_name: Name of output vcf file
  :param chrom: Chromosome number
  :param ml: master list of variants
  :param samples: list of samples
  :param sample_names: list of strings representing sample names
  :param bgzip_and_index: If True (default)
  """
  with open(vcf_name, 'w') as fp:
    # Write header
    sample_header = '\t'.join(sample_names)
    fp.write(
      "##fileformat=VCFv4.1\n"
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{:s}\n".format(sample_header)
    )
    # Write lines
    wr = fp.write
    viter = get_variant_rows(ml, samples)
    # TODO


def sort_and_index_bam(bamfile):
  """Do the filename gymnastics required to end up with a sorted, indexed, bam file."""
  # samtools sort adds a '.bam' to the end of the file name.
  os.rename(bamfile, 'temp.bam')
  pysam.sort('temp.bam', os.path.splitext(bamfile)[0])
  pysam.index(bamfile)
  os.remove('temp.bam')


def vcf2chrom(vcf_rdr):
  """Given a vcf reader corresponding to one chromosome, read in the data descriptions into our format. The result is
  sorted if the vcf file is sorted.
  """
  chrom = []  # deque()
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

    append(new_variation(start, stop, ref, alt, het))

  return chrom


def parse_vcf(vcf_rdr, chrom_list):
  """Given a vcf reader load in all the chromosomes for all the samples."""
  g1 = {}
  for chrom in chrom_list:
    try:
      g1[chrom] = vcf2chrom(vcf_rdr.fetch(chrom, start=0))
    except (ValueError, KeyError):  # New version of pyvcf changed the error
      g1[chrom] = []  #deque()
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
  for chrom, variants in g1.iteritems():
    ch = str(chrom)
    for var in variants:
      # In the VCF file no REF or ALT is indicated by a .
      ref = var.vd.REF if var.vd.REF != '' else '.'
      alt = var.vd.ALT if var.vd.ALT != '' else '.'
      #  CHROM    POS   ID   REF   ALT   QUAL FILTER INFO FORMAT tsample
      wr(ch + "\t" + str(var.vd.POS) + "\t.\t" + ref + "\t" + alt + "\t100\tPASS\t.\tGT\t" + GT[var.het] + "\n")


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


def compress_and_index_vcf(in_vcf_name, out_vcf_name):
  """Given an uncompressed, but sorted, vcf, compress and index it."""
  #bgzip -c sorted.vcf > sorted.vcf.gz
  #tabix sorted.vcf.gz
  logger.debug('Compressing and indexing {:s} to {:s}'.format(in_vcf_name, out_vcf_name))
  pysam.tabix_compress(in_vcf_name, out_vcf_name, force=True)
  pysam.tabix_index(out_vcf_name, force=True, preset='vcf')


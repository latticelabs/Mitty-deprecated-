"""Some utilities related to file formats, compressing and indexing."""
from os.path import splitext
import os
import cPickle
import sqlite3 as sq
from contextlib import contextmanager

import pysam

#from mitty.lib.variation import HET_01, HET_10, HOMOZYGOUS, ABSENT, new_variation, GT
from mitty.lib.variation import VariationData as VD

import logging
logger = logging.getLogger(__name__)


@contextmanager
def open_db(db_name='population.sqlite3'):
  conn = sq.connect(db_name)
  conn.text_factory = str
  c = conn.cursor()
  yield conn, c
  conn.close()


def db(db_name='population.sqlite3'):
  conn = sq.connect(db_name)
  conn.text_factory = str
  return conn


def save_variant_master_list(chrom_name, ml, conn=None):
  """Save the master list as a chromosome in a sqlite3 database
  :param chrom_name: name of the chromosome. The table will be named
  :param ml: dictionary representing the variant master list
  :param conn: The connection object.
  """
  table_name = 'chrom_' + chrom_name
  conn.execute("DROP TABLE IF EXISTS {:s}".format(table_name))
  conn.execute("CREATE TABLE {:s} (vid INTEGER PRIMARY KEY, idx INTEGER, pos INTEGER, ref TEXT, alt TEXT)".format(table_name))
  for k, v in ml.iteritems():
    conn.execute("INSERT INTO {:s}(idx, pos, ref, alt) VALUES (?, ?, ?, ?)".format(table_name), (k, v.POS, v.REF, v.ALT))
  conn.commit()


def load_variant_master_list(chrom_name, conn):
  """Save the master list as a python pickle file
  :param chrom_name: name of the chromosome. The table will be named
  :param conn: The connection object.
  :returns ml: dictionary representing the variant master list
  """
  table_name = 'chrom_' + chrom_name
  ml = {}
  for row in conn.execute("SELECT * FROM {:s} ORDER BY idx".format(table_name)):
    ml[row[1]] = VD(row[2], row[2] + len(row[3]), row[3], row[4])
  return ml


def sort_and_index_bam(bamfile):
  """Do the filename gymnastics required to end up with a sorted, indexed, bam file."""
  # samtools sort adds a '.bam' to the end of the file name.
  os.rename(bamfile, 'temp.bam')
  pysam.sort('temp.bam', os.path.splitext(bamfile)[0])
  pysam.index(bamfile)
  os.remove('temp.bam')


def vcf2chrom(vcf_rdr):
  """Given a vcf reader corresponding to one chromosome, read in the variant descriptions into our format. The result is
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


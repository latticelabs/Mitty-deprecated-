import vcf
import io
from mitty.lib.variation import *
from . import *  # To get definitions from the setup script
from nose.tools import assert_equal


def vcf2chrom_test():

  #
  #  INS     DEL--     SNP      DEL---     INV-----
  #  0    1  2 3 4  5  6    7   8 9 10  11 12 13 14
  #  ------  -----     ------   ------     --------

  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t1\t.\tC\tCAA\t100\tPASS\t.\tGT\t0/1\n"
    "1\t3\t.\tCAG\tC\t100\tPASS\t.\tGT\t1/0\n"
    "1\t7\t.\tG\tT\t100\tPASS\t.\tGT\t0/1\n"
    "1\t9\t.\tGTT\t.\t100\tPASS\t.\tGT\t1/0\n"
    "1\t13\t.\tGTT\tTTG\t100\tPASS\t.\tGT\t1/1\n"
  )

  correct_chrom = [
      new_variation(1, 2, 'C', 'CAA', HET_01),
      new_variation(3, 6, 'CAG', 'C', HET_10),
      new_variation(7, 8, 'G', 'T', HET_01),
      new_variation(9, 12, 'GTT', '', HET_10),
      new_variation(13, 16, 'GTT', 'TTG', HOMOZYGOUS)
  ]
  chrom = vcf2chrom(vcf.Reader(fsock=io.BytesIO(vcf_str)))
  assert chrom == correct_chrom, chrom


def vcf2chrom_test2():
  """Read VCF file with no sample/genotype information."""

  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "1\t1\t.\tC\tCAA\t100\tPASS\t.\n"
    "1\t3\t.\tCAG\tC\t100\tPASS\t.\n"
    "1\t7\t.\tG\tT\t100\tPASS\t.\n"
    "1\t9\t.\tGTT\t.\t100\tPASS\t.\n"
    "1\t13\t.\tGTT\tTTG\t100\tPASS\t.\n"
  )

  correct_chrom = [
      new_variation(1, 2, 'C', 'CAA', HOMOZYGOUS),
      new_variation(3, 6, 'CAG', 'C', HOMOZYGOUS),
      new_variation(7, 8, 'G', 'T', HOMOZYGOUS),
      new_variation(9, 12, 'GTT', '', HOMOZYGOUS),
      new_variation(13, 16, 'GTT', 'TTG', HOMOZYGOUS)
  ]
  chrom = vcf2chrom(vcf.Reader(fsock=io.BytesIO(vcf_str)))
  assert chrom == correct_chrom, chrom


def vcf_empty_chromosome_test():
  """Ask for chromosome that does not exist in VCF file"""
  g1 = parse_vcf(vcf.Reader(filename=small_vcf_name + '.gz'), [1, 2, 3])
  assert g1[3] == []


def vcf_round_trip_test():
  """VCF round trip (load and then save)."""
  import tempfile
  g1 = parse_vcf(vcf.Reader(filename=small_vcf_name + '.gz'), [1, 2])
  temp_vcf_fp, temp_vcf_name = tempfile.mkstemp(suffix='.vcf')  # No .gz extension on purpose
  vcf_save_gz(g1, temp_vcf_name)
  g1_load = parse_vcf(vcf.Reader(filename=temp_vcf_name + '.gz'), [1, 2])
  assert_equal(g1, g1_load)
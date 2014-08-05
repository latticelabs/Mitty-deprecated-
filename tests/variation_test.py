from mitty.variation import *
from . import *


def vcf2chrom_test():
  import vcf, io

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
      Variation(1, 2, 'C', 'CAA', HET2),
      Variation(3, 6, 'CAG', 'C', HET1),
      Variation(7, 8, 'G', 'T', HET2),
      Variation(9, 12, 'GTT', '', HET1),
      Variation(13, 16, 'GTT', 'TTG', HOMOZYGOUS)
  ]
  chrom = vcf2chrom(vcf.Reader(fsock=io.BytesIO(vcf_str)))
  assert chrom == correct_chrom, chrom


def vcf2chrom_test2():
  """Read VCF file with no sample/genotype information."""
  import vcf, io

  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "1\t1\t.\tC\tCAA\t100\tPASS\t.\n"
    "1\t3\t.\tCAG\tC\t100\tPASS\t.\n"
    "1\t7\t.\tG\tT\t100\tPASS\t.\n"
    "1\t9\t.\tGTT\t.\t100\tPASS\t.\n"
    "1\t13\t.\tGTT\tTTG\t100\tPASS\t.\n"
  )

  correct_chrom = [
      Variation(1, 2, 'C', 'CAA', HOMOZYGOUS),
      Variation(3, 6, 'CAG', 'C', HOMOZYGOUS),
      Variation(7, 8, 'G', 'T', HOMOZYGOUS),
      Variation(9, 12, 'GTT', '', HOMOZYGOUS),
      Variation(13, 16, 'GTT', 'TTG', HOMOZYGOUS)
  ]
  chrom = vcf2chrom(vcf.Reader(fsock=io.BytesIO(vcf_str)))
  assert chrom == correct_chrom, chrom


def vcf_round_trip_test():
  """VCF round trip (load and then save)."""
  import io
  original = open(small_vcf_name, 'r').readlines()
  g1 = parse_vcf(vcf.Reader(filename=small_vcf_name + '.gz'), [1, 2])
  fp = io.BytesIO()
  vcf_save(g1, fp)
  fp.seek(0)
  ln = fp.readlines()[1:]  # The first line is the VCF version which we ignore
  assert original == ln, [original, ln]

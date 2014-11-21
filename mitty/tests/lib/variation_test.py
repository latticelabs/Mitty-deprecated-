import vcf
import io
from mitty.lib.variation import *
from .. import *  # To get definitions from the setup script
from nose.tools import assert_sequence_equal, assert_sequence_equal


def vcf2chrom_test():
  """Loading from VCF formatted string."""
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
  assert_sequence_equal(chrom, correct_chrom)


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
  assert_sequence_equal(chrom, correct_chrom)


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
  assert_sequence_equal(g1, g1_load)


def merge_test1():
  """Merge variants, non overlapping, existing first (ED)."""
  c1 = [v10] = [new_variation(1, 4, 'CAA', 'C', HOMOZYGOUS)]
  c2 = [v20] = [new_variation(10, 13, 'CAA', 'C', HOMOZYGOUS)]
  c3 = merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10, v20])


def merge_test2():
  """Merge variants, non overlapping, denovo first (DE)."""
  c1 = [v10] = [new_variation(10, 13, 'CAA', 'C', HOMOZYGOUS)]
  c2 = [v20] = [new_variation(1, 4, 'CAA', 'C', HOMOZYGOUS)]
  c3 = merge_variants(c1, c2)
  assert_sequence_equal(c3, [v20, v10])


def merge_test3():
  """Merge variants, non overlapping (EDED)"""
  c1 = [v10, v11] = [new_variation(1, 4, 'CAA', 'C', HOMOZYGOUS), new_variation(13, 16, 'CTT', 'C', HOMOZYGOUS)]
  c2 = [v20, v21] = [new_variation(8, 11, 'CCC', 'C', HOMOZYGOUS), new_variation(20, 23, 'CAA', 'C', HOMOZYGOUS)]
  c3 = merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10, v20, v11, v21])


def merge_test4a():
  """Merge variants, full overlapping (E-D-)"""
  c1 = [v10] = [new_variation(2, 5, 'CAA', 'C', HOMOZYGOUS)]
  c2 = [v20] = [new_variation(2, 5, 'CAA', 'T', HOMOZYGOUS)]
  c3 = merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10])


def merge_test4b():
  """Merge variants, full overlapping, SNP (D-E-)."""
  c1 = [v10] = [new_variation(2, 3, 'C', 'G', HET_10)]
  c2 = [v20] = [new_variation(2, 3, 'C', 'T', HET_10)]
  c3 = merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10])


# Important test - killed a nasty logic bug
def merge_test4c():
  """Merge variants, full overlapping, with non-colliding preceder"""
  c1 = [v10, v11] = [new_variation(1, 2, 'G', 'C', HET_01),
                     new_variation(2, 3, 'A', 'C', HET_10)]
  c2 = [v20] = [new_variation(2, 3, 'A', 'T', HET_10)]
  c3 = merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10, v11])


def merge_test4():
  """Merge variants, overlapping (E-D). D will collide and will be rejected"""
  c1 = [v10] = [new_variation(1, 4, 'CAA', 'C', HOMOZYGOUS)]
  c2 = [v20] = [new_variation(2, 5, 'CCC', 'C', HOMOZYGOUS)]
  c3 = merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10])


def merge_test5():
  """Merge variants, overlapping (D-E). D will collide and will be rejected"""
  c2 = [v20] = [new_variation(1, 4, 'CAA', 'C', HOMOZYGOUS)]
  c1 = [v10] = [new_variation(2, 5, 'CCC', 'C', HOMOZYGOUS)]
  c3 = merge_variants(c1, c2)
  assert_sequence_equal(c3, [v10])


def merge_test6():
  """Merge variants, overlapping heterozygous. Overlapping but with mixed zygosity"""
  c1 = [v10, v11, v12, v13] = [new_variation(1, 4, 'CAA', 'C', HOMOZYGOUS),
                               new_variation(13, 16, 'CAA', 'C', HET_10),
                               new_variation(20, 23, 'CAA', 'C', HET_10),
                               new_variation(26, 29, 'CAA', 'C', HOMOZYGOUS)]
  c2 = [v20, v21, v22, v23] = [new_variation(5, 8, 'CTT', 'C', HOMOZYGOUS),  # Will collide out
                               new_variation(13, 16, 'CTT', 'C', HET_01),  # Will pass
                               new_variation(20, 23, 'CTT', 'C', HOMOZYGOUS),  # Will collide out
                               new_variation(26, 29, 'CTT', 'C', HET_01)]  # Will collide out
  c3 = merge_variants(c1, c2)
  # assert_sequence_equal(c2, [new_variation(1, 4, 'CAA', 'C', HOMOZYGOUS),
  #                   new_variation(13, 16, 'CTT', 'C', HET_10),
  #                   new_variation(13, 16, 'CAA', 'C', HET_01),
  #                   new_variation(20, 23, 'CAA', 'C', HET_10),
  #                   new_variation(26, 29, 'CGG', 'C', HOMOZYGOUS)])
  assert_sequence_equal(deque([v10, v11, v21, v12, v13]), c3)


def merge_test7():
  """Skip bug, discovered during simulations"""
  cv1 = new_variation(3, 7, 'ACTG', 'A', HOMOZYGOUS)
  cv2 = new_variation(10, 11, 'C', 'A', HOMOZYGOUS)
  c1 = [cv1, cv2]

  dv1 = new_variation(1, 5, 'ACTG', 'A', HOMOZYGOUS)
  dv2 = new_variation(7, 8, 'C', 'T', HOMOZYGOUS)  # Collides with previous, should be discarded
  dv3 = new_variation(15, 16, 'A', 'T', HOMOZYGOUS)  # Too close together, second one should collide out
  dv4 = new_variation(16, 16, 'C', 'T', HOMOZYGOUS)
  dnv = [dv1, dv2, dv3, dv4]

  c2 = merge_variants(c1, dnv)
  assert_sequence_equal(c2, [cv1, cv2, dv3])

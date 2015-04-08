import tempfile
import os

import vcf

import mitty.lib.io as mio
import mitty.lib.variants as vr
from mitty.tests import *  # To get definitions from the setup script
from nose.tools import assert_raises
#from nose.tools import assert_list_equal, assert_dict_equal


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
  chrom = [(0, 0), (2, 1), (3,2)]

  with mio.vcf_for_writing(temp_name, ['a']) as fp:
    mio.write_chromosomes_to_vcf(fp, seq_id='chr2', chrom_list=[chrom], master_list=ml)

  rdr = vcf.Reader(filename=temp_name)
  v = rdr.next()
  assert v.POS == 2
  assert v.REF == 'A'
  assert v.ALT[0].sequence == 'AA'
  assert v.genotype('a').gt_type == 1

  v = rdr.next()
  v = rdr.next()
  assert v.POS == 31
  assert v.REF == 'GAAAA'
  assert v.ALT[0].sequence == 'G'
  assert v.genotype('a').gt_type == 2



# def master_list_round_trip_test():
#   """Save a master list and then load it back."""
#   l = {}
#   sample = [vr.cgt(l, 1, 4, 'CAA', 'C', vr.HOMOZYGOUS),
#             vr.cgt(l, 13, 16, 'CAA', 'C', vr.HET_10),
#             vr.cgt(l, 20, 23, 'CAA', 'C', vr.HET_10),
#             vr.cgt(l, 26, 29, 'CAA', 'C', vr.HOMOZYGOUS)]
#   temp_fp, temp_name = tempfile.mkstemp(suffix='.sqlite3')
#   os.close(temp_fp)
#   conn = mio.db(temp_name)
#   mio.save_variant_master_list('chr1', l, conn)
#   rt_l = mio.load_variant_master_list('chr1', conn)
#   assert_dict_equal(l, rt_l, rt_l)
#   os.remove(temp_name)
#
#
# def sample_list_round_trip_test():
#   """Save a sample and then load it back."""
#   l = {}
#   sample = [vr.cgt(l, 1, 4, 'CAA', 'C', vr.HOMOZYGOUS),
#             vr.cgt(l, 13, 16, 'CAA', 'C', vr.HET_10),
#             vr.cgt(l, 20, 23, 'CAA', 'C', vr.HET_10),
#             vr.cgt(l, 26, 29, 'CAA', 'C', vr.HOMOZYGOUS)]
#   temp_fp, temp_name = tempfile.mkstemp(suffix='.sqlite3')
#   os.close(temp_fp)
#   conn = mio.db(temp_name)
#   mio.save_sample('s1', 'chr1', sample, conn)
#   rt_sample = mio.load_sample('s1', 'chr1', conn)
#   #mio.save_variant_master_list('chr1', l, conn)
#   #rt_l = mio.load_variant_master_list('chr1', conn)
#   assert_list_equal(sample, rt_sample, rt_sample)
#   os.remove(temp_name)


# def vcf2chrom_test():
#   """Loading from VCF formatted string."""
#   #
#   #  INS     DEL--     SNP      DEL---     INV-----
#   #  0    1  2 3 4  5  6    7   8 9 10  11 12 13 14
#   #  ------  -----     ------   ------     --------
#
#   vcf_str = (
#     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
#     "1\t1\t.\tC\tCAA\t100\tPASS\t.\tGT\t0/1\n"
#     "1\t3\t.\tCAG\tC\t100\tPASS\t.\tGT\t1/0\n"
#     "1\t7\t.\tG\tT\t100\tPASS\t.\tGT\t0/1\n"
#     "1\t9\t.\tGTT\t.\t100\tPASS\t.\tGT\t1/0\n"
#     "1\t13\t.\tGTT\tTTG\t100\tPASS\t.\tGT\t1/1\n"
#   )
#
#   correct_chrom = [
#       new_variation(1, 2, 'C', 'CAA', HET_01),
#       new_variation(3, 6, 'CAG', 'C', HET_10),
#       new_variation(7, 8, 'G', 'T', HET_01),
#       new_variation(9, 12, 'GTT', '', HET_10),
#       new_variation(13, 16, 'GTT', 'TTG', HOM)
#   ]
#   chrom = vcf2chrom(vcf.Reader(fsock=io.BytesIO(vcf_str)))
#   assert_sequence_equal(chrom, correct_chrom)
#
#
# def vcf2chrom_test2():
#   """Read VCF file with no sample/genotype information."""
#
#   vcf_str = (
#     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
#     "1\t1\t.\tC\tCAA\t100\tPASS\t.\n"
#     "1\t3\t.\tCAG\tC\t100\tPASS\t.\n"
#     "1\t7\t.\tG\tT\t100\tPASS\t.\n"
#     "1\t9\t.\tGTT\t.\t100\tPASS\t.\n"
#     "1\t13\t.\tGTT\tTTG\t100\tPASS\t.\n"
#   )
#
#   correct_chrom = [
#       new_variation(1, 2, 'C', 'CAA', HOM),
#       new_variation(3, 6, 'CAG', 'C', HOM),
#       new_variation(7, 8, 'G', 'T', HOM),
#       new_variation(9, 12, 'GTT', '', HOM),
#       new_variation(13, 16, 'GTT', 'TTG', HOM)
#   ]
#   chrom = vcf2chrom(vcf.Reader(fsock=io.BytesIO(vcf_str)))
#   assert_sequence_equal(chrom, correct_chrom)
#
#
# def vcf_empty_chromosome_test():
#   """Ask for chromosome that does not exist in VCF file"""
#   g1 = parse_vcf(vcf.Reader(filename=small_vcf_name + '.gz'), [1, 2, 3])
#   assert g1[3] == []
#
#
# def vcf_round_trip_test():
#   """VCF round trip (load and then save)."""
#   import tempfile
#   g1 = parse_vcf(vcf.Reader(filename=small_vcf_name + '.gz'), [1, 2])
#   temp_vcf_fp, temp_vcf_name = tempfile.mkstemp(suffix='.vcf')  # No .gz extension on purpose
#   vcf_save_gz(g1, temp_vcf_name)
#   g1_load = parse_vcf(vcf.Reader(filename=temp_vcf_name + '.gz'), [1, 2])
#   assert_sequence_equal(g1, g1_load)

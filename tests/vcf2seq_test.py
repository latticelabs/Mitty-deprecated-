import vcf
import pysam
import io
import numpy.testing
from . import *
from mitty import vcf2seq
from nose.tools import assert_equals, ok_
from shutil import rmtree

def assembly_check(ref_seq, vcf_str,
                   var_seq_1, pos_1, var_seq_2, pos_2,
                   v_coords, v_codes):
  """Utility function to test if two sequence assemblies and their pos vectors are the same."""
  def check(tup1, tup2):
    copy_1, p_1, v_pos_1, v_cod_1 = tup1
    copy_2, p_2, v_pos_2, v_cod_2 = tup2
    for n in [0,1]:
      assert_equals(len(copy_1[n]), len(copy_2[n]))
      assert_equals(len(p_1[n]), len(p_2[n]))

      for x, y in zip(copy_1[n], copy_2[n]):
        assert_equals(x, y)

      numpy.testing.assert_array_equal(p_1[n], p_2[n])

    numpy.testing.assert_array_equal(v_pos_1, v_pos_2, err_msg=str(v_pos_1) + '|' + str(v_pos_2))
    numpy.testing.assert_array_equal(v_cod_1, v_cod_2)

  a2 = ([var_seq_1, var_seq_2], [pos_1, pos_2], v_coords, v_codes)
  a1 = vcf2seq.assemble_sequences(ref_seq, vcf.Reader(fsock=io.BytesIO(vcf_str)))
  check(a1, a2)


def test_assemble_sequences_homo_snp():
  """Testing homozygous SNP assembly"""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\tT\t100\tPASS\t.\tGT\t1/1\n"
  )
  var_seq_1 = 'ATTG'
  pos_1 = [1, 2, 3, 4, 5]
  var_seq_2 = 'ATTG'
  pos_2 = [1, 2, 3, 4, 5]
  v_coords = [[[1, 2], [1, 2]]]
  v_codes = [[1, 0]]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2, v_coords, v_codes)


def test_assemble_sequences_hetero_snp():
  """Testing heterozygous SNP assembly"""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\tT\t100\tPASS\t.\tGT\t1/0\n"
  )
  var_seq_1 = 'ATTG'
  pos_1 = [1, 2, 3, 4, 5]
  var_seq_2 = 'ACTG'
  pos_2 = [1, 2, 3, 4, 5]
  v_coords = [[[1, 2], [1, 1]]]
  v_codes = [[1, 1]]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2, v_coords, v_codes)


def test_assemble_sequences_hetero_ins():
  """Testing heterozygous insertion assembly"""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\tCC\t100\tPASS\t.\tGT\t1/0\n"
  )
  var_seq_1 = 'ACCTG'
  pos_1 = [1, 2, 3, 3, 4, 5]
  var_seq_2 = 'ACTG'
  pos_2 = [1, 2, 3, 4, 5]
  v_coords = [[[1, 3], [1, 1]]]
  v_codes = [[2, 1]]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2, v_coords, v_codes)


def test_assemble_sequences_hetero_del():
  """Testing heterozygous deletion assembly C -> ."""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\t.\t100\tPASS\t.\tGT\t0/1\n"
  )
  var_seq_1 = 'ACTG'
  pos_1 = [1, 2, 3, 4, 5]
  var_seq_2 = 'ATG'
  pos_2 = [1, 3, 4, 5]
  v_coords = [[[1, 1], [1, 1]]]
  v_codes = [[2, 2]]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2, v_coords, v_codes)


def test_assemble_sequences_hetero_del2():
  """Testing heterozygous deletion assembly CT -> C"""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tCT\tC\t100\tPASS\t.\tGT\t0/1\n"
  )
  var_seq_1 = 'ACTG'
  pos_1 = [1, 2, 3, 4, 5]
  var_seq_2 = 'ACG'
  pos_2 = [1, 2, 4, 5]
  v_coords = [[[1, 1], [1, 2]]]
  v_codes = [[2, 2]]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2, v_coords, v_codes)


def test_assemble_sequences_hetero_same_locus_del():
  """Testing repeated VCF entries (same location, same type - SNP - different chromosome copies)"""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\t.\t100\tPASS\t.\tGT\t0/1\n"
    "1\t2\t.\tC\t.\t100\tPASS\t.\tGT\t1/0\n"
  )
  var_seq_1 = 'ATG'
  pos_1 = [1, 3, 4, 5]
  var_seq_2 = 'ATG'
  pos_2 = [1, 3, 4, 5]
  v_coords = [[[1, 1], [1, 1]], [[1, 1], [1, 1]]]
  v_codes = [[2, 2], [2, 1]]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2, v_coords, v_codes)


def test_assemble_sequences_hetero_same_locus_ins_del():
  """Testing repeated VCF entries (same location, same type - SNP, different formats - different chromosome copies)"""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\t.\t100\tPASS\t.\tGT\t0/1\n"
    "1\t2\t.\tC\tCT\t100\tPASS\t.\tGT\t1/0\n"
  )
  var_seq_1 = 'ACTTG'
  pos_1 = [1, 2, 3, 3, 4, 5]
  var_seq_2 = 'ATG'
  pos_2 = [1, 3, 4, 5]
  v_coords = [[[1, 1], [1, 1]], [[1, 3], [1, 1]]]
  v_codes = [[2, 2], [2, 1]]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2, v_coords, v_codes)


def test_assemble_sequences_hetero_same_locus_ins_del2():
  """Testing repeated VCF entries (same location, one del, one ins, different chromosome copies)"""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tCTG\t.\t100\tPASS\t.\tGT\t0/1\n"
    "1\t2\t.\tC\tCT\t100\tPASS\t.\tGT\t1/0\n"
  )
  var_seq_1 = 'ACTTG'
  pos_1 = [1, 2, 3, 3, 4, 5]
  var_seq_2 = 'A'
  pos_2 = [1, 5]
  v_coords = [[[1, 1], [1, 1]], [[1, 3], [1, 1]]]
  v_codes = [[2, 2], [2, 1]]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2, v_coords, v_codes)


def test_assemble_sequences_hetero_same_locus_inv_del():
  """Testing repeated VCF entries (same location, deletion and inversion, different chromosome copies)"""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tCT\tGA\t100\tPASS\t.\tGT\t0/1\n"
    "1\t2\t.\tC\t.\t100\tPASS\t.\tGT\t1/0\n"
  )
  var_seq_1 = 'ATG'
  pos_1 = [1, 3, 4, 5]
  var_seq_2 = 'AGAG'
  pos_2 = [1, 2, 3, 4, 5]

  v_coords = [[[1, 1], [1, 3]], [[1, 1], [3, 3]]]
  v_codes = [[2, 2], [2, 1]]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2, v_coords, v_codes)


def test_assemble_sequences_hetero_same_locus_longer_ins():
  """Testing multiple repeated VCF entries (heterozygous insertions)"""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\tCAA\t100\tPASS\t.\tGT\t0/1\n"
    "1\t2\t.\tC\tCGG\t100\tPASS\t.\tGT\t1/0\n"
    "1\t4\t.\tG\tGAA\t100\tPASS\t.\tGT\t0/1\n"
    "1\t4\t.\tG\tGTT\t100\tPASS\t.\tGT\t1/0\n"
  )
  var_seq_1 = 'ACGGTGTT'
  pos_1 = [1, 2, 3, 3, 3, 4, 5, 5, 5]
  var_seq_2 = 'ACAATGAA'
  pos_2 = [1, 2, 3, 3, 3, 4, 5, 5, 5]
  v_coords = [[[1, 1], [1, 4]], [[1, 4], [4, 4]], [[5, 5], [5, 8]], [[5, 8], [8 ,8]]]
  v_codes = [[2, 2], [2, 1], [2, 2], [2, 1]]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2, v_coords, v_codes)


def test_script():
  """vcf2seq  command line program"""

  ok_(os.path.exists(wg_name),
      msg='No whole genome file ({:s}). This should be created by package test setup in tests/__init__.py'.format(wg_name))

  tempdir = tempfile.mkdtemp()
  vcf_name = os.path.join(tempdir, 'vcf2seq_test.vcf')
  var_name = os.path.join(tempdir, 'vcf2seq_test_var.h5')

  args = {
    '--ref': wg_name,  # Our package wide setup as generated this
    '--vcf': vcf_name + '.gz',
    '--var': var_name
  }

  # Insert in copy 2, delete in copy 1
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t71\t.\tG\tGCAT\t100\tPASS\t.\tGT\t0/1\n"
    "1\t71\t.\tGGCG\tG\t100\tPASS\t.\tGT\t1/0\n"
  )

  with open(vcf_name, 'w') as fp:
    fp.write(vcf_str)

  pysam.tabix_compress(vcf_name, vcf_name + '.gz', force=True)
  pysam.tabix_index(vcf_name + '.gz', force=True, preset='vcf')

  vcf2seq.main(args)

  ok_(os.path.exists(var_name), msg='Output file was not created')

  #Now check some details of the file that convinces us vcf2seq is working as expected
  with h5py.File(var_name, 'r') as h5_fp:
    assert not h5_fp['sequence/1/1'].attrs['reference']
    assert not h5_fp['sequence/1/2'].attrs['reference']
    assert h5_fp['sequence/2/2'].attrs['reference']
    assert_equals(h5_fp['sequence/1/2'][70:74].tostring(), 'GCAT')
    assert_equals(h5_fp['sequence/1/1'][70:74].tostring(), 'GGCG')

  rmtree(tempdir)  # Be neat
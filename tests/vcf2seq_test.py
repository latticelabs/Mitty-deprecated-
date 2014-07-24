import vcf
import io
from mitty import vcf2seq
from nose.tools import assert_equals


def assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, variant_coords_1, var_seq_2, pos_2, variant_coords_2):
  """Utility function to test if two sequence assemblies and their pos vectors are the same."""
  def check(tup1, tup2):
    copy_1, p_1, vc_dict_1 = tup1
    copy_2, p_2, vc_dict_2 = tup2
    for n in [0,1]:
      assert_equals(len(copy_1[n]), len(copy_2[n]))
      assert_equals(len(p_1[n]), len(p_2[n]))

      for x, y in zip(copy_1[n], copy_2[n]):
        assert_equals(x, y)

      for x, y in zip(p_1[n], p_2[n]):
        assert_equals(x, y)

      assert_equals(vc_dict_1[n], vc_dict_2[n])

  a2 = ([var_seq_1, var_seq_2], [pos_1, pos_2], [variant_coords_1, variant_coords_2])
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
  vc_1 = {'ins': []}
  var_seq_2 = 'ATTG'
  pos_2 = [1, 2, 3, 4, 5]
  vc_2 = {'ins': []}

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, vc_1, var_seq_2, pos_2, vc_2)


def test_assemble_sequences_hetero_snp():
  """Testing heterozygous SNP assembly"""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\tT\t100\tPASS\t.\tGT\t1/0\n"
  )
  var_seq_1 = 'ATTG'
  pos_1 = [1, 2, 3, 4, 5]
  vc_1 = {'ins': []}
  var_seq_2 = 'ACTG'
  pos_2 = [1, 2, 3, 4, 5]
  vc_2 = {'ins': []}

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, vc_1, var_seq_2, pos_2, vc_2)


def test_assemble_sequences_hetero_ins():
  """Testing heterozygous insertion assembly"""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\tCC\t100\tPASS\t.\tGT\t1/0\n"
  )
  var_seq_1 = 'ACCTG'
  pos_1 = [1, 2, 3, 3, 4, 5]
  vc_1 = {'ins': [[1, 3]]}
  var_seq_2 = 'ACTG'
  pos_2 = [1, 2, 3, 4, 5]
  vc_2 = {'ins': []}

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, vc_1, var_seq_2, pos_2, vc_2)


def test_assemble_sequences_hetero_del():
  """Testing heterozygous deletion assembly C -> ."""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\t.\t100\tPASS\t.\tGT\t0/1\n"
  )
  var_seq_1 = 'ACTG'
  pos_1 = [1, 2, 3, 4, 5]
  vc_1 = {'ins': []}
  var_seq_2 = 'ATG'
  pos_2 = [1, 3, 4, 5]
  vc_2 = {'ins': []}

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, vc_1, var_seq_2, pos_2, vc_2)


def test_assemble_sequences_hetero_del2():
  """Testing heterozygous deletion assembly CT -> C"""
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tCT\tC\t100\tPASS\t.\tGT\t0/1\n"
  )
  var_seq_1 = 'ACTG'
  pos_1 = [1, 2, 3, 4, 5]
  vc_1 = {'ins': []}
  var_seq_2 = 'ACG'
  pos_2 = [1, 2, 4, 5]
  vc_2 = {'ins': []}

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, vc_1, var_seq_2, pos_2, vc_2)


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
  vc_1 = {'ins': []}
  var_seq_2 = 'ATG'
  pos_2 = [1, 3, 4, 5]
  vc_2 = {'ins': []}

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, vc_1, var_seq_2, pos_2, vc_2)


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
  vc_1 = {'ins': [[1,3]]}
  var_seq_2 = 'ATG'
  pos_2 = [1, 3, 4, 5]
  vc_2 = {'ins': []}

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, vc_1, var_seq_2, pos_2, vc_2)


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
  vc_1 = {'ins': [[1,3]]}
  var_seq_2 = 'A'
  pos_2 = [1, 5]
  vc_2 = {'ins': []}

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, vc_1, var_seq_2, pos_2, vc_2)


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
  vc_1 = {'ins': []}
  var_seq_2 = 'AGAG'
  pos_2 = [1, 2, 3, 4, 5]
  vc_2 = {'ins': [[1, 3]]}  # This is interesting, because vcf flags this as an insertion

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, vc_1, var_seq_2, pos_2, vc_2)


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
  vc_1 = {'ins': [[1,4], [5, 8]]}
  var_seq_2 = 'ACAATGAA'
  pos_2 = [1, 2, 3, 3, 3, 4, 5, 5, 5]
  vc_2 = {'ins': [[1, 4], [5, 8]]}

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, vc_1, var_seq_2, pos_2, vc_2)


def test_assemble_sequences_misc():
  """Testing a more comprehensive VCF file"""
  ref_seq = 'ACTGACTGACTG'
  vcf_str = (
    "##fileformat=VCFv4.1\n"
    "##fileDate={:s}\n"
    "##source=mutate.py {:s} ({:s})\n"
    "##reference={:s}\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t1\t.\tA\tT\t100\tPASS\t.\tGT\t0/1\n"
    "1\t2\t.\tC\tT\t100\tPASS\t.\tGT\t1/0\n"
    "1\t3\t.\tT\tG\t100\tPASS\t.\tGT\t1/1\n"
    "1\t4\t.\tG\t.\t100\tPASS\t.\tGT\t0/1\n"
    "1\t5\t.\tAC\tA\t100\tPASS\t.\tGT\t1/0\n"
    "1\t6\t.\tCT\tC\t100\tPASS\t.\tGT\t0/1\n"
    "1\t8\t.\tG\tGC\t100\tPASS\t.\tGT\t0/1\n"
    "1\t8\t.\tG\tGG\t100\tPASS\t.\tGT\t1/0\n"
  )
  var_seq_1 = 'ATGGATGGACTG'
  pos_1 = [1,2,3,4,5,7,8,9,9,10,11,12, 13]
  vc_1 = {'ins': [[6, 8]]}
  var_seq_2 = 'TCGACGCACTG'
  pos_2 = [1,2,3,5,6,8,9,9,10,11,12, 13]
  vc_2 = {'ins': [[5, 7]]}

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, vc_1, var_seq_2, pos_2, vc_2)



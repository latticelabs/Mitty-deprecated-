import vcf
import io
from .. import vcf2seq
from nose.tools import assert_equals


def assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2):
  def check(a1, a2):
    copy_a, pos_a = a1
    copy_b, pos_b = a2
    for n in [0,1]:
      assert_equals(len(copy_a[n]), len(copy_b[n]))
      assert_equals(len(pos_a[n]), len(pos_b[n]))

    for n in [0,1]:
      for x, y in zip(copy_a[n], copy_b[n]):
        assert_equals(x, y)

      for x, y in zip(pos_a[n], pos_b[n]):
        assert_equals(x, y)

  a2 = ([var_seq_1, var_seq_2], [pos_1, pos_2])
  a1 = vcf2seq.assemble_sequences(ref_seq, vcf.Reader(fsock=io.BytesIO(vcf_str)))
  check(a1, a2)


def test_assemble_sequences_1():
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\tT\t100\tPASS\t.\tGT\t1/1\n"
  )
  var_seq_1 = 'ATTG'
  pos_1 = [1, 2, 3, 4, 5]
  var_seq_2 = 'ATTG'
  pos_2 = [1, 2, 3, 4, 5]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2)


def test_assemble_sequences_2():
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\tT\t100\tPASS\t.\tGT\t1/0\n"
  )
  var_seq_1 = 'ATTG'
  pos_1 = [1, 2, 3, 4, 5]
  var_seq_2 = 'ACTG'
  pos_2 = [1, 2, 3, 4, 5]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2)


def test_assemble_sequences_3():
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\tCC\t100\tPASS\t.\tGT\t1/0\n"
  )
  var_seq_1 = 'ACCTG'
  pos_1 = [1, 2, 3, 3, 4, 5]
  var_seq_2 = 'ACTG'
  pos_2 = [1, 2, 3, 4, 5]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2)


def test_assemble_sequences_4():
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tC\t.\t100\tPASS\t.\tGT\t0/1\n"
  )
  var_seq_1 = 'ACTG'
  pos_1 = [1, 2, 3, 4, 5]
  var_seq_2 = 'ATG'
  pos_2 = [1, 3, 4, 5]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2)


def test_assemble_sequences_5():
  ref_seq = 'ACTG'
  vcf_str = (
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    "1\t2\t.\tCT\tC\t100\tPASS\t.\tGT\t0/1\n"
  )
  var_seq_1 = 'ACTG'
  pos_1 = [1, 2, 3, 4, 5]
  var_seq_2 = 'ACG'
  pos_2 = [1, 2, 4, 5]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2)


def test_assemble_sequences_6():
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

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2)


def test_assemble_sequences_7():
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

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2)


def test_assemble_sequences_8():
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

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2)


def test_assemble_sequences_9():
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

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2)


def test_assemble_sequences_10():
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
  var_seq_2 = 'TCGACGCACTG'
  pos_2 = [1,2,3,5,6,8,9,9,10,11,12, 13]

  assembly_check(ref_seq, vcf_str, var_seq_1, pos_1, var_seq_2, pos_2)

"""Test cases for BADBAM intersection and difference. It is suffiicent to generate a test case for one chromosome
as the badbams code operates on reads on a chromosome by chromosome basis"""
import tempfile
import os

import pysam
import numpy as np

from mitty.benchmarking.badbams import intersect_and_diff_of_bams
from mitty.lib.mio import sort_and_index_bam


def write_badbam(bam_fp, chrom, cpy_array, ro_array, pos_array, pos_off=1, index=[], rl=100, tl=250):
  """Given a list of pos values, write this out to a BADBAM file. All errors are chrom + pos
  (with both being incremented by one). We create reads misaligned in both files but misaligned
  differently by changing the pos_off value from its default.

  TAG TYPE VALUE
  Zc  i    0 - read comes from chrom copy 0, 1 - read comes from chrom copy 1
  ZE  i    Read stop (Read start is in POS)
  Ze  i    Mate stop (Mate start is available from other BAM info)
  Xf  i    0 - incorrectly mapped, 1 - correctly mapped, 2 - unmapped
  YR  i    0 - chrom was wrong, 1 - chrom was correct
  YP  i    0 - pos was wrong, 1 - pos was correct
  YC  i    0 - CIGAR was wrong, 1 - CIGAR was correct
  XR  i    Aligned chromosome
  XP  i    Aligned pos
  XC  Z    Aligned CIGAR

  @read_serial|chrom|copy|ro|pos|rlen|cigar|ro|pos|rlen|cigar

  """
  cigar = '{:d}M'.format(rl)
  #for n, cpy, ro, pos in enumerate(zip(cpy_array, ro_array, pos_array)):
  for n in index:
    cpy, ro, pos = cpy_array[n], ro_array[n], pos_array[n]
    p1, p2 = pos, pos + tl
    new_read = pysam.AlignedSegment()
    new_read.reference_id = chrom
    new_read.pos = pos
    new_read.qname = '{read_serial}|{chrom}|{copy}|{ro1}|{pos1}|{rlen}|{cigar}|{ro2}|{pos2}|{rlen}|{cigar}'.\
      format(read_serial=n, chrom=chrom, copy=cpy, ro1=ro, pos1=p1, rlen=rl, cigar=cigar, ro2=1-ro, pos2=p2)

    # Needs to be consistent with __extended_bam_tags_info__
    new_read.set_tags([('Zc', cpy, 'i'),
                       ('ZE', pos + rl, 'i'),
                       ('Ze', pos + tl + rl, 'i'),
                       ('Xf', 0, 'i'),
                       ('YR', 0, 'i'),
                       ('YP', 0, 'i'),
                       ('YC', 0, 'i'),
                       ('XR', cpy + 1, 'i'),
                       ('XP', pos + pos_off, 'i'),
                       ('XC', cigar, 'Z')])
    bam_fp.write(new_read)


def create_two_badbams(bam1, bam2):
  """Create two BADBAM files with known characteristics.

  Make a simple list of reads. Assign 33% to be misaligned in both BADBAMs, 33% to be misaligned in one and 33% to be
  misaligned in the other.

  Of the 33% misaligned in both, the first 50% are identical and last 50% are differently misaligned
  """
  hdr = {
    'HD': {'SO': 'coordinate', 'VN': '1.3'},
    'SQ': [{'SN': 'chrom1', 'LN': 100}, {'SN': 'chrom2', 'LN': 100}]
  }
  rng = np.random.RandomState(seed=123456)
  n_reads, max_pos, chrom = 99, 50, 1
  pos = np.sort(np.random.randint(max_pos, size=n_reads))
  cpy = np.sort(np.random.randint(2, size=n_reads))
  ro = np.sort(np.random.randint(2, size=n_reads))  # Strand of first read

  chrom = 0
  shf_idx = rng.permutation(n_reads)
  t0 = int(n_reads / 6.0)
  t1 = int(n_reads / 3.0)
  tt = [0, t0, t1, int(n_reads / 1.5), n_reads]
  idx_both_identical, idx_both_diff, idx_1, idx_2 = shf_idx[tt[0]:tt[1]], shf_idx[tt[1]:tt[2]], shf_idx[tt[2]:tt[3]], shf_idx[tt[3]:tt[4]]

  for n, fname in enumerate([bam1, bam2]):
    fp = pysam.AlignmentFile(fname, 'wb', header=hdr)
    # Identically misaligned reads
    write_badbam(fp, chrom, cpy, ro, pos, index=idx_both_identical, rl=100, tl=250)
    # Differently misaligned reads in both
    write_badbam(fp, chrom, cpy, ro, pos, index=idx_both_diff, pos_off=n + 1, rl=100, tl=250)
    # Misaligned reads present in one or the other
    diff_idx = idx_1 if n == 0 else idx_2
    write_badbam(fp, chrom, cpy, ro, pos, index=diff_idx, rl=100, tl=250)
    fp.close()
    # Sort and index
    sort_and_index_bam(fname)

  return chrom, cpy, ro, pos, idx_both_identical, idx_both_diff, idx_1, idx_2, hdr


def analyse_two_badbams(bam1, bam2, hdr):
  _, intersect_bam_identical_errors = tempfile.mkstemp(suffix='.bam')
  _, intersect_bam_diff_errors_1 = tempfile.mkstemp(suffix='.bam')
  _, intersect_bam_diff_errors_2 = tempfile.mkstemp(suffix='.bam')
  _, diff_bam1 = tempfile.mkstemp(suffix='.bam')
  _, diff_bam2 = tempfile.mkstemp(suffix='.bam')

  intersect_bam_identical_errors_fp = pysam.AlignmentFile(intersect_bam_identical_errors, 'wb', header=hdr)
  intersect_bam_diff_errors_1_fp = pysam.AlignmentFile(intersect_bam_diff_errors_1, 'wb', header=hdr)
  intersect_bam_diff_errors_2_fp = pysam.AlignmentFile(intersect_bam_diff_errors_2,'wb', header=hdr)
  diff_bam1_fp = pysam.AlignmentFile(diff_bam1, 'wb', header=hdr)
  diff_bam2_fp = pysam.AlignmentFile(diff_bam2, 'wb', header=hdr)

  bam_fp1 = pysam.AlignmentFile(bam1, 'rb')
  bam_fp2 = pysam.AlignmentFile(bam2, 'rb')

  for _ in intersect_and_diff_of_bams(
    bam_fp1=bam_fp1, bam_fp2=bam_fp2,
    intersect_bam_identical_errors=intersect_bam_identical_errors_fp,
    intersect_bam_diff_errors_1=intersect_bam_diff_errors_1_fp,
    intersect_bam_diff_errors_2=intersect_bam_diff_errors_2_fp,
    diff_bam1=diff_bam1_fp, diff_bam2=diff_bam2_fp):
    pass

  return intersect_bam_identical_errors, intersect_bam_diff_errors_1, intersect_bam_diff_errors_2, diff_bam1, diff_bam2


def assert_bam_is_correct(bam, index):
  """We simply load the list of read serials, and verify they are the same as that required from the index"""
  read_index = [int(r.qname.split('|')[0]) for r in pysam.AlignmentFile(bam, 'rb')]
  assert len(read_index) > 0, 'No reads in this set'
  assert set(read_index) == set(index), read_index


def test_badbam():
  """Test badbam operation"""
  _, bam1 = tempfile.mkstemp(dir='./', suffix='.bam')
  _, bam2 = tempfile.mkstemp(dir='./', suffix='.bam')
  chrom, cpy, ro, pos, idx_both_identical, idx_both_diff, idx_1, idx_2, hdr = create_two_badbams(bam1, bam2)
  intersect_bam_identical_errors, intersect_bam_diff_errors_1, intersect_bam_diff_errors_2, diff_bam1, diff_bam2 = analyse_two_badbams(bam1, bam2, hdr)

  # Now read in the files and make sure they contain what they should contain
  assert_bam_is_correct(intersect_bam_identical_errors, idx_both_identical)
  assert_bam_is_correct(intersect_bam_diff_errors_1, idx_both_diff)
  assert_bam_is_correct(intersect_bam_diff_errors_2, idx_both_diff)
  assert_bam_is_correct(diff_bam1, idx_1)
  assert_bam_is_correct(diff_bam2, idx_2)

  # Cleanup
  os.remove(bam1)
  os.remove(bam1 + '.bai')
  os.remove(bam2)
  os.remove(bam2 + '.bai')
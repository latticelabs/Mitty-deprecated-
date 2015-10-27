"""CLI and library for taking in a pair of BADBAM files and looking at intersect and difference. Saves summary as .json
file and writes out intersect and difference files as BADBAMs and FASTQs. The FASTQs provide small datasets for
investigating a problem and the subset BADBAMs are useful for analysing what kind of errors were made"""
import time
import io

import click
import pysam

from mitty.version import __version__

import logging
logger = logging.getLogger(__name__)


def bams_have_same_reference(bam_fp1, bam_fp2):
  h1 = bam_fp1.header['SQ']
  h2 = bam_fp2.header['SQ']
  if len(h1) != len(h2): return False
  for s1, s2 in zip(h1, h2):
    if s1['LN'] != s2['LN']: return False
  return True


def read_errors_are_identical(r1, r2):
  """r1 and r2 have been identified as the same read. Are their errors identical?

  Xf  i    0 - incorrectly mapped, 1 - correctly mapped, 2 - unmapped
  YR  i    0 - chrom was wrong, 1 - chrom was correct
  YP  i    0 - pos was wrong, 1 - pos was correct
  YC  i    0 - CIGAR was wrong, 1 - CIGAR was correct
  XR  i    Aligned chromosome
  XP  i    Aligned pos
  XC  Z    Aligned CIGAR
  """
  for tag in ['XC', 'XP', 'XR']:
    if r1.get_tag(tag) != r2.get_tag(tag): return False
  return True


def intersect_and_diff_of_bams(bam_fp1, bam_fp2,
                               intersect_bam_identical_errors,
                               intersect_bam_diff_errors_1,
                               intersect_bam_diff_errors_2,
                               diff_bam1, diff_bam2,
                               progress_update_interval=1000):
  """Given two position sorted BADBAMs compute the intersection and difference of the two

  :param bam_fp1: first bam file
  :param bam_fp2: second bam file
  :param intersect_bam_identical_errors: same reads with identical errors
  :param intersect_bam_diff_errors_1: same reads with diff errors - read from bam1 goes here
  :param intersect_bam_diff_errors_2: same reads with diff errors - read from bam2 goes here
  :param diff_bam1: reads only wrong in bam1
  :param diff_bam2: reads only wrong in bam2
  :param progress_update_interval: every N reads pass a message back (via yield) to indicate how many reads we've done

  Algorithm:
  Advance along the two BAMs.
  If pos1 < pos2 then r1 is part of 1 - 2. Advance fp1
  If pos2 < pos2 then r2 is part of 2 - 1. Advance fp2
  If pos1 = pos2,
      if qnames match then r1,r2 are part of 1 ^ 2,
      else r1 is part of 1 -2, r2 is part of 2 -1
      advance both
  """
  cntr = progress_update_interval
  for seq_id in [sn['SN'] for sn in bam_fp1.header['SQ']]:
    itr1, itr2 = bam_fp1.fetch(reference=seq_id), bam_fp2.fetch(reference=seq_id)
    r1, r2 = next(itr1, None), next(itr2, None)
    while r1 is not None and r2 is not None:
      if r1.pos < r2.pos:
        diff_bam1.write(r1)
        r1 = next(itr1, None)
      elif r2.pos < r1.pos:
        diff_bam2.write(r2)
        r2 = next(itr2, None)
      elif r1.qname != r2.qname:
        diff_bam1.write(r1)
        diff_bam2.write(r2)
        r1, r2 = next(itr1, None), next(itr2, None)
      else:
        # If the correct positions are the same and the qnames are the same, they can not be mates so
        # we don't have to test for that
        if read_errors_are_identical(r1, r2):
          intersect_bam_identical_errors.write(r1)
        else:
          intersect_bam_diff_errors_1.write(r1)
          intersect_bam_diff_errors_2.write(r2)
        r1, r2 = next(itr1, None), next(itr2, None)
      cntr -= 1
      if cntr == 0:
        yield progress_update_interval
        cntr = progress_update_interval

    # Pickup the slack

    cntr = progress_update_interval
    while r1 is not None:
      diff_bam1.write(r1)
      r1 = next(itr1, None)
      cntr -= 1
      if cntr == 0:
        yield progress_update_interval
        cntr = progress_update_interval

    cntr = progress_update_interval
    while r2 is not None:
      diff_bam2.write(r2)
      r2 = next(itr2, None)
      cntr -= 1
      if cntr == 0:
        yield progress_update_interval
        cntr = progress_update_interval

    yield progress_update_interval  # This is a white lie (To save time we aren't counting the actual number of reads)


@click.command()
@click.version_option()
@click.argument('inbam1', type=click.Path(exists=True))
@click.argument('inbam2', type=click.Path(exists=True))
@click.option('--intersect-bam-identical-errors', type=click.Path(), help='same reads with identical errors')
@click.option('--intersect_bam_diff_errors_1', type=click.Path(), help='same reads with diff errors - read from bam1 goes here')
@click.option('--intersect_bam_diff_errors_2', type=click.Path(), help='same reads with diff errors - read from bam2 goes here')
@click.option('--diff_bam1', type=click.Path(), help='reads only wrong in bam1')
@click.option('--diff_bam2', type=click.Path(), help='reads only wrong in bam2')
@click.option('-v', count=True, help='Verbosity level')
@click.option('-p', is_flag=True, help='Show progress bar')
def cli(inbam1, inbam2,
        intersect_bam_identical_errors,
        intersect_bam_diff_errors_1,
        intersect_bam_diff_errors_2,
        diff_bam1, diff_bam2, v, p):
  """Take a pair of BADBAM files and look at intersect and difference"""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  bam_fp1 = pysam.AlignmentFile(inbam1, 'rb')
  bam_fp2 = pysam.AlignmentFile(inbam2, 'rb')

  if not bams_have_same_reference(bam_fp1, bam_fp2):
    logger.error('Reference file in BAMs are not identical. Exiting')
    exit(0)

  new_header = bam_fp1.header
  new_header['PG'].append({
    'CL': 'badbams',
    'ID': 'mitty-badbams',
    'PN': 'badbam',
    'VN': __version__,
    'PP': new_header['PG'][-1]['ID']
  })

  intersect_bam_identical_errors_fp = pysam.AlignmentFile(intersect_bam_identical_errors or 'intersect_identical.bam',
                                                          'wb', header=new_header)
  intersect_bam_diff_errors_1_fp = pysam.AlignmentFile(intersect_bam_diff_errors_1 or 'intersect_diff_errors1.bam',
                                                       'wb', header=new_header)
  intersect_bam_diff_errors_2_fp = pysam.AlignmentFile(intersect_bam_diff_errors_2 or 'intersect_diff_errors2.bam',
                                                       'wb', header=new_header)
  diff_bam1_fp = pysam.AlignmentFile(diff_bam1 or 'diff_bam1.bam', 'wb', header=new_header)
  diff_bam2_fp = pysam.AlignmentFile(diff_bam2 or 'diff_bam2.bam', 'wb', header=new_header)

  cnt = 0
  t0 = time.time()
  total_read_count = max(bam_fp1.mapped + bam_fp1.unmapped, bam_fp2.mapped + bam_fp2.unmapped)  # Sadly, this is only approximate
  progress_bar_update_interval = int(0.01 * total_read_count)
  cnt = 0
  with click.progressbar(length=total_read_count, label='Processing BAMs',
                         file=None if p else io.BytesIO()) as bar:
    for _ in intersect_and_diff_of_bams(
      bam_fp1=bam_fp1, bam_fp2=bam_fp2,
      intersect_bam_identical_errors=intersect_bam_identical_errors_fp,
      intersect_bam_diff_errors_1=intersect_bam_diff_errors_1_fp,
      intersect_bam_diff_errors_2=intersect_bam_diff_errors_2_fp,
      diff_bam1=diff_bam1_fp, diff_bam2=diff_bam2_fp,
      progress_update_interval=progress_bar_update_interval):
      bar.update(progress_bar_update_interval)
      cnt += progress_bar_update_interval
  t1 = time.time()
  logger.debug('Analyzed {:d} reads in BAM in {:2.2f}s'.format(cnt, t1 - t0))
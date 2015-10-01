"""Go through a BAM file made from alignments on a FASTQ and do the following:

1. BADBAM: Create a bam file containing only the misaligned reads. The read (CHROM, POS, CIGAR) are set to be the correct alignment
and the original alignment (CHROM, POS, CIGAR) is encoded in extended tags (see below).
The reason for setting the read CHROM, POS and CIGAR attributes to the correct one is to be able to merge several of
these BAM files. Having the reads in the correct (and therefore fixed) order enables us to merge sort them

2. PERBAM: Create a BAM file containing all reads with the correct (CHROM, POS, CIGAR) set and the aligned (CHROM,POS,CIGAR)
stored using extended tags (same as above, details below).

If we are asked for just a read analysis then this file omits read data (such as sequence and quality scores) so that
the resultant BAM is smaller and more manageable.

If we are asked for a perfect BAM then we also write the full read data into this file and can use this as a perfect
input for variant callers and so on.

Extended tags

TAG TYPE VALUE
Zc  A    0 - read comes from chrom copy 0, 1 - read comes from chrom copy 1
Xf  A    0 - incorrectly mapped, 1 - correctly mapped, 2 - unmapped
YR  A    0 - chrom was wrong, 1 - chrom was correct
YP  A    0 - pos was wrong, 1 - pos was correct
YC  A    0 - CIGAR was wrong, 1 - CIGAR was correct
XR  i    Aligned chromosome
XP  i    Aligned pos
XC  Z    Aligned CIGAR"""
import os
import time
import io

import pysam
import click

from mitty.version import __version__
import mitty.benchmarking.creed as creed
import mitty.lib.io as mio  # For the bam sort and index function
from mitty.lib import DNA_complement
from string import translate

import logging
logger = logging.getLogger(__name__)


def process_file(bam_in_fp, bad_bam_fp=None, per_bam_fp=None, full_perfect_bam=False, window=0, extended=False,
                 flag_cigar_errors_as_misalignments=False,
                 progress_bar_update_interval=100):
  """Main processing function that goes through the bam file, analyzing read alignment and writing out

  :param bam_in_fp:  Pointer to original BAM
  :param bad_bam_fp: Pointer to BADBAM being created
  :param per_bam_fp: Pointer to PERBAM being created
  :param full_perfect_bam: If True, write the read seq and qual scores too, to make a complete perfect bam
  :param window:     Tolerance window for deciding if read is correctly aligned
  :param extended:   If True write out new style CIGARs (With '=' and 'X')
  :param flag_cigar_errors_as_misalignments: Set to True if we want CIGAR errors to count as misalignments
  :param progress_bar_update_interval: how many reads to process before yielding (to update progress bar as needed)
  :return: number of reads processed
  """
  n0 = progress_bar_update_interval
  analyze_read = creed.analyze_read
  for cnt, read in enumerate(bam_in_fp):
    read_serial, chrom, cpy, ro, pos, cigar, chrom_c, pos_c, cigar_c, unmapped = analyze_read(read, window, extended)
    #if read_serial != 8200: continue
    if read_serial is None: continue  # Something wrong with this read.
    read_is_misaligned = not (chrom_c and pos_c and (cigar_c or (not flag_cigar_errors_as_misalignments)))
    if read_is_misaligned or full_perfect_bam:  # Need all the read info, incl seq and quality
      new_read = read
    else:  # Need only some tags and pos, chrom info
      new_read = pysam.AlignedSegment()

    # Zc  A    0 - read comes from chrom copy 0, 1 - read comes from chrom copy 1
    # Xf  A    0 - incorrectly mapped, 1 - correctly mapped, 2 - unmapped
    # YR  A    0 - chrom was wrong, 1 - chrom was correct
    # YP  A    0 - pos was wrong, 1 - pos was correct
    # YC  A    0 - CIGAR was wrong, 1 - CIGAR was correct
    # XR  i    Aligned chromosome
    # XP  i    Aligned pos
    # XC  Z    Aligned CIGAR
    new_read.set_tags([('Zc', str(cpy), 'A'),
                       ('Xf', '2' if unmapped else str(chrom_c and pos_c), 'A'),
                       ('YR', str(chrom_c), 'A'),
                       ('YP', str(pos_c), 'A'),
                       ('YC', str(cigar_c), 'A'),
                       ('XR', read.reference_id, 'i'),
                       ('XP', read.pos, 'i'),
                       ('XC', read.cigarstring or '', 'Z')])

    if read_is_misaligned:  # We need to check if the read was reverse complemented when it should have been and vv
      if new_read.is_reverse != ro:  # The complement is not consistent
        qual = new_read.qual[::-1]
        new_read.seq = translate(new_read.seq, DNA_complement)[::-1]
        new_read.qual = qual

    new_read.is_reverse = ro
    new_read.mate_is_reverse = 1 - ro
    new_read.mate_is_unmapped = False  # Gotta check this - what if mate is deep in an insert?
    new_read.reference_id = chrom - 1
    new_read.pos = pos
    new_read.cigarstring = cigar  # What if this is deep in an insert?

    if read_is_misaligned:
      bad_bam_fp.write(new_read)

    per_bam_fp.write(new_read)

    n0 -= 1
    if n0 == 0:
      yield cnt
      n0 = progress_bar_update_interval


def analyze(in_bam_fname, bad_bam_fname, per_bam_fname, flag_cigar_errors, perfect_bam, window, x, p):
  bam_in_fp = pysam.AlignmentFile(in_bam_fname, 'rb')

  def true2str(v): return 'true' if v else 'false'

  new_header = bam_in_fp.header
  new_header['PG'].append({
    'CL': 'perfectbam ....',
    'ID': 'mitty-perfectbam',
    'PN': 'perfectbam',
    'VN': __version__,
    'PP': new_header['PG'][-1]['ID'],
    'DS': 'window={:d}, cigar errors result in misalignments={:s}, extended_cigar={:s}'.
      format(window, true2str(flag_cigar_errors), true2str(x))
  })

  bad_bam_fp = pysam.AlignmentFile(bad_bam_fname, 'wb', header=new_header)
  per_bam_fp = pysam.AlignmentFile(per_bam_fname, 'wb', header=new_header)

  cnt = 0
  t0 = time.time()
  total_read_count = bam_in_fp.mapped + bam_in_fp.unmapped  # Sadly, this is only approximate
  progress_bar_update_interval = int(0.01 * total_read_count)
  with click.progressbar(length=total_read_count, label='Processing BAM',
                         file=None if p else io.BytesIO()) as bar:
    for cnt in process_file(bam_in_fp=bam_in_fp, bad_bam_fp=bad_bam_fp, per_bam_fp=per_bam_fp,
                            full_perfect_bam=perfect_bam, window=window,
                            flag_cigar_errors_as_misalignments=flag_cigar_errors, extended=x,
                            progress_bar_update_interval=progress_bar_update_interval):
      bar.update(progress_bar_update_interval)
  t1 = time.time()
  logger.debug('Analyzed {:d} reads in BAM in {:2.2f}s'.format(cnt, t1 - t0))


def sort_and_index(bad_bam_fname, per_bam_fname):
  t0 = time.time()
  mio.sort_and_index_bam(bad_bam_fname)
  t1 = time.time()
  logger.debug('Sort and indexed bad BAM in {:2.2f}s'.format(t1 - t0))

  t0 = time.time()
  mio.sort_and_index_bam(per_bam_fname)
  t1 = time.time()
  logger.debug('Sort and indexed perfect BAM in {:2.2f}s'.format(t1 - t0))


@click.command()
@click.version_option()
@click.argument('inbam', type=click.Path(exists=True))
@click.option('--cigar-errors-are-misalignments', is_flag=True, help='CIGAR errors result in reads being classified as misaligned')
@click.option('--perfect-bam', is_flag=True, help='Write out perfect BAM')
@click.option('--window', help='Size of tolerance window', default=0, type=int)
@click.option('-x', is_flag=True, help='Use extended CIGAR ("X"s and "="s) rather than traditional CIGAR (just "M"s)')
@click.option('-v', count=True, help='Verbosity level')
@click.option('-p', is_flag=True, help='Show progress bar')
def cli(inbam, cigar_errors_are_misalignments, perfect_bam, window, x, v, p):
  """Analyse BAMs produced from Mitty generated fastqs for alignment accuracy."""
  level = logging.DEBUG if v > 0 else logging.WARNING
  logging.basicConfig(level=level)

  bad_bam_fname = os.path.splitext(inbam)[0] + '_bad.bam'
  per_bam_fname = os.path.splitext(inbam)[0] + '_per.bam'

  analyze(inbam, bad_bam_fname, per_bam_fname, cigar_errors_are_misalignments, perfect_bam, window, x, p)
  sort_and_index(bad_bam_fname, per_bam_fname)


if __name__ == "__main__":
  cli()
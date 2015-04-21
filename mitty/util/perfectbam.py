#!python
__cmd__ = """Commandline::

  Usage:
    perfectbam <inbam> [<out_prefix>] [--window=WN] [-x] [-v] [-p]

  Options:
    <inbam>        Input bam file name of reads
    <out_prefix>   Written files will have this prefix. If absent, the prefix will be the same as the bam file
    --window=WN    Size of tolerance window [default: 0]
    -x             Use extended CIGAR ('X's and '='s) rather than traditional CIGAR (just 'M's)
    -v             Dump detailed logger messages
    -p             Show progress bar
"""
__param__ = """Given a bam file containing simulated reads aligned by a tool:
  1. Produce a new bam that re-aligns all reads so that their alignment is perfect
  2. Produce a csv file containing data about the misaligned and unmapped reads with the following columns::
       qname              -  qname of the read
       error_type         -  type of error 3 bit number  bit 0=chrom, 1=pos, 2=cigar
       correct_chrom      -  correct chromosome number of read
       correct_pos        -  correct position of read
       correct_cigar      -  correct cigar of read
       aligned_chrom      -  actual aligned chromosome
       aligned_pos        -  actual aligned pos
       aligned_cigar      -  actual aligned cigar
       mapping_qual       -  mapping quality
       mate_is_unmapped   -  is mate unmapped
       seq                -  actual sequence string
  3. Produce a .json file with various useful summary statistics"""
__doc__ = __cmd__ + __param__

import os
import json
import time

import pysam
from collections import Counter
import docopt

import mitty.lib.io as mio  # For the bam sort and index function
from mitty.lib.reads import old_style_cigar
from mitty.lib import progress_bar

import logging
logger = logging.getLogger(__name__)


def csv_header():
  header = [
    'qname',
    'error_type',
    'correct_chrom',
    'correct_pos',
    'correct_cigar',
    'aligned_chrom',
    'aligned_pos',
    'aligned_cigar',
    'mapping_qual',
    'mate_is_unmapped',
    'seq',
  ]
  return header


def main(bam_in_fp, bam_out_fp, csv_fp, json_fp, window, extended=False, progress_bar_func=None):
  """Main processing function that goes through the bam file, analyzing read alignment and writing out

  :param bam_in_fp:
  :param bam_out_fp:
  :param csv_fp:
  :param json_fp:
  :param window:
  :return: number of reads processed
  """
  csv_fp.write('\t'.join(csv_header()) + '\n')

  total_read_count = float(bam_in_fp.mapped + bam_in_fp.unmapped)
  f0 = 0
  total_reads_cntr, incorrectly_aligned_reads_cntr = Counter(), Counter()

  if progress_bar_func is not None: progress_bar_func('Processing BAM ', 0, 80)

  for n, read in enumerate(bam_in_fp):
    # qname = 'r{:d}|{:d}|{:d}|{:d}|{:s}|{:d}|{:s}'
    if read.is_paired:
      if read.is_read1:
        rd_ser, chrom, cpy, ro, pos, cigar, _, _, _ = read.qname.split('|')
      else:
        rd_ser, chrom, cpy, _, _, _, ro, pos, cigar = read.qname.split('|')
    else:
      rd_ser, chrom, cpy, ro, pos, cigar = read.qname.split('|')
    ro, chrom, pos = int(ro), int(chrom), int(pos)
    if not extended:
      cigar = old_style_cigar(cigar)

    total_reads_cntr[chrom] += 1

    error_type = 0x0
    if read.reference_id != chrom - 1:
      error_type |= 0x1
    if not (-window <= read.pos - pos <= window):
      error_type |= 0x2
    if read.cigarstring != cigar:
      error_type |= 0x4
    if error_type != 0:
      incorrectly_aligned_reads_cntr[chrom] += 1
      csv_fp.write(', '.join(map(str, [read.qname, error_type, chrom, pos, cigar,
                                       read.reference_id + 1, read.pos, read.cigarstring,
                                       read.mapq, read.mate_is_unmapped, read.query_sequence])) + '\n')

    # Now write out the perfect alignment
    read.is_reverse = 1 - ro
    read.mate_is_reverse = ro
    read.mate_is_unmapped = False  # Gotta check this - what if mate is deep in an insert?
    read.reference_id = chrom - 1
    read.pos = pos
    read.cigarstring = cigar  # What if this is deep in an insert?
    bam_out_fp.write(read)

    if progress_bar_func is not None:
      f = n / total_read_count
      if f - f0 >= 0.01:
        progress_bar_func('Processing BAM ', f, 80)
        f0 = f
  if progress_bar_func is not None: print('\n')

  json.dump({'sequence_header': bam_in_fp.header['SQ'],
             'read_counts': {str(k): v for k,v in total_reads_cntr.iteritems()},
             'incorrectly_aligned_read_counts': {str(k): v for k, v in incorrectly_aligned_reads_cntr.iteritems()}},
            json_fp, indent=2)

  return int(total_read_count)


def cli():
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__cmd__, ['-h'])
  else:
    args = docopt.docopt(__cmd__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  fname_prefix = args['<out_prefix>'] or os.path.splitext(args['<inbam>'])[0]
  csv_fname = fname_prefix + '_misaligned.csv'
  summary_fname = fname_prefix + '_summary.json'
  perfect_bam_fname = fname_prefix + '_perfect.bam'

  with pysam.AlignmentFile(args['<inbam>'], 'rb') as bam_in_fp, \
      pysam.AlignmentFile(perfect_bam_fname, 'wb', template=bam_in_fp) as bam_out_fp, \
      open(csv_fname, 'w') as csv_out_fp, open(summary_fname, 'w') as json_out_fp:
    t0 = time.time()
    read_count = main(bam_in_fp=bam_in_fp, bam_out_fp=bam_out_fp, csv_fp=csv_out_fp, json_fp=json_out_fp,
                      window=int(args['--window']), extended=bool(args['-x']), progress_bar_func=progress_bar if args['-p'] else None)
    t1 = time.time()
    logger.debug('Analyzed {:d} reads in BAM in {:2.2f}s'.format(read_count, t1 - t0))
  t0 = time.time()
  mio.sort_and_index_bam(perfect_bam_fname)
  t1 = time.time()
  logger.debug('Sort and indexed perfect BAM in {:2.2f}s'.format(t1 - t0))


if __name__ == "__main__":
  cli()
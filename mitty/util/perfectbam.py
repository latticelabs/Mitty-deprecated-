#!python
__cmd__ = """Commandline::

  Usage:
    perfectbam  --inbam=INBAM  [--window=WN] [-v]

  Options:
    --inbam=INBAM           Input bam file name of reads
    --window=WN             Size of tolerance window [default: 0]
    -v                      Dump detailed logger messages
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
import sys

import pysam
from collections import Counter
import docopt

import mitty.lib.io as mio  # For the bam sort and index function

import logging
logger = logging.getLogger(__name__)


# def analyze_bam(bam_fp, window=0, block_size=100000):
#   """Given a read iterator from a BAM file, process each read and return a list of lists suitable for writing to csv
#   file."""
#   total_reads_cntr, bad_reads_cntr = Counter(), Counter()
#   keep_reading = True
#   while keep_reading:
#     misaligned_reads = [None] * block_size
#     current_count = block_size
#     while current_count > 0:
#       try:
#         read = bam_fp.next()
#         correct_chrom_no, _, correct_pos, correct_cigar = interpret_read_qname(read.qname, read.is_read2)
#         total_reads_cntr[correct_chrom_no] += 1
#         if correct_chrom_no == read.tid + 1 and abs(correct_pos - (read.pos + 1)) <= window:  # BAM is zero indexed
#           # This read is correctly aligned
#           continue
#         misaligned_reads[block_size - current_count] = \
#           [read.qname, correct_chrom_no, correct_pos, read.tid + 1, read.pos + 1, read.mapq, read.mate_is_unmapped, read.seq]
#         bad_reads_cntr[correct_chrom_no] += 1
#         current_count -= 1
#       except StopIteration:
#         keep_reading = False
#         break
#     if current_count > 0:
#       del misaligned_reads[block_size - current_count:]
#
#     yield misaligned_reads, total_reads_cntr, bad_reads_cntr


# def write_csv_header(csv_fp):
#   csv_fp.write('qname, correct_chrom, correct_pos, aligned_chrom, aligned_pos, mapping_qual, mate_is_unmapped, seq\n')
#
#
# def write_csv(csv_fp, misaligned_reads):
#   csv_fp.writelines((', '.join(map(str, line)) + '\n' for line in misaligned_reads))


def main(bam_in_fp, bam_out_fp, csv_fp, json_fp, window):
  """Main processing function that goes through the bam file, analyzing read alignment and writing out

  :param bam_in_fp:
  :param bam_out_fp:
  :param csv_fp:
  :param json_fp:
  :param window:
  :return:
  """
  total_read_count = float(bam_in_fp.count())
  bam_in_fp.reset()  # This is bad pysam design. The count() action iterates through the file!
  f0 = 0
  total_reads_cntr, incorrectly_aligned_reads_cntr = Counter(), Counter()
  for n, read in enumerate(bam_in_fp):
    # qname = 'r{:d}|{:d}|{:d}|{:d}|{:s}|{:d}|{:s}'
    if read.is_paired:
      if read.is_read1:
        rd_ser, chrom, cpy, pos, cigar, _, _ = read.qname.split('|')
      else:
        rd_ser, chrom, cpy, _, _, pos, cigar = read.qname.split('|')
    else:
      rd_ser, chrom, cpy, pos, cigar = read.qname.split('|')
    chrom, pos = int(chrom), int(pos)
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
    read.reference_id = chrom - 1
    read.pos = pos
    read.cigarstring = cigar  # What if this is deep in an insert?
    read.mate_is_unmapped = False  # Gotta check this - what if mate is deep in an insert?
    bam_out_fp.write(read)

    f = n / total_read_count
    if f - f0 >= 0.01:
      progress_bar('Processing BAM', f, 80)
      f0 = f
  print('\n')

  json.dump({"read_counts": {str(k): v for k,v in total_reads_cntr.iteritems()},
             "incorrectly_aligned_read_counts": {str(k): v for k, v in incorrectly_aligned_reads_cntr.iteritems()}},
            json_fp, indent=2)


def progress_bar(title, f, cols):
  """Draw a nifty progress bar.
  '\r' trick from http://stackoverflow.com/questions/15685063/print-a-progress-bar-processing-in-python

  :param title: leading text to print
  :param f:     fraction completed
  :param cols:  how many columns wide should the bar be
  """
  x = int(f * cols + 0.5)
  sys.stdout.write('\r' + title + '[' + '.' * x + ' ' * (cols - x) + ']')
  sys.stdout.flush()


def cli():
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__cmd__, ['-h'])
  else:
    args = docopt.docopt(__cmd__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  fname_prefix = os.path.splitext(args['--inbam'])[0]
  csv_fname = fname_prefix + '_misaligned.csv'
  summary_fname = fname_prefix + '_summary.json'
  perfect_bam_fname = fname_prefix + '_perfect.bam'

  with pysam.AlignmentFile(args['--inbam'], 'rb') as bam_in_fp, \
      pysam.AlignmentFile(perfect_bam_fname, 'wb', template=bam_in_fp) as bam_out_fp, \
      open(csv_fname, 'w') as csv_out_fp, open(summary_fname, 'w') as json_out_fp:
    main(bam_in_fp=bam_in_fp, bam_out_fp=bam_out_fp, csv_fp=csv_out_fp, json_fp=json_out_fp,
         window=int(args['--window']))
  mio.sort_and_index_bam(perfect_bam_fname)


if __name__ == "__main__":
  cli()
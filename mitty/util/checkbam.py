#!python
"""Given a bam file containing simulated reads aligned by a tool produce diagnostic files for analyzing the aligner
performance:
  1. Produce a csv file containing data about the misaligned and unmapped reads with the following columns::
       qname, correct_chrom, correct_pos, aligned_chrom, aligned_pos, mapping_qual, mate_is_unmapped, seq
  2. Produce a .json file with various useful summary statistics

Commandline::

  Usage:
    checkbam  --inbam=INBAM  --out_prefix=OUT_PREF  [--window=WN] [--block_size=BS] [-v]

  Options:
    --inbam=INBAM           Input bam file name of reads
    --out_prefix=OUT_PREF   Prefix of output files
    --window=WN             Size of tolerance window [default: 0]
    --block_size=BS         Number of reads to process at a time [default: 1000000]
    -v                      Dump detailed logger messages
"""
__version__ = '1.0.0'

import pysam
import json
from collections import Counter
import docopt
from mitty.vcf2reads import interpret_read_qname
import logging
logger = logging.getLogger(__name__)


def analyze_bam(bam_fp, window=0, block_size=100000):
  """Given a read iterator from a BAM file, process each read and return a list of lists suitable for writing to csv
  file."""
  total_reads_cntr, bad_reads_cntr = Counter(), Counter()
  keep_reading = True
  while keep_reading:
    misaligned_reads = [None] * block_size
    current_count = block_size
    while current_count > 0:
      try:
        read = bam_fp.next()
        correct_chrom_no, _, correct_pos, correct_cigar = interpret_read_qname(read.qname, read.is_read2)
        total_reads_cntr[correct_chrom_no] += 1
        if correct_chrom_no == read.tid + 1 and abs(correct_pos - (read.pos + 1)) <= window:  # BAM is zero indexed
          # This read is correctly aligned
          continue
        misaligned_reads[block_size - current_count] = \
          [read.qname, correct_chrom_no, correct_pos, read.tid + 1, read.pos + 1, read.mapq, read.mate_is_unmapped, read.seq]
        bad_reads_cntr[correct_chrom_no] += 1
        current_count -= 1
      except StopIteration:
        keep_reading = False
        break
    if current_count > 0:
      del misaligned_reads[block_size - current_count:]

    yield misaligned_reads, total_reads_cntr, bad_reads_cntr


def write_csv_header(csv_fp):
  csv_fp.write('qname, correct_chrom, correct_pos, aligned_chrom, aligned_pos, mapping_qual, mate_is_unmapped, seq\n')


def write_csv(csv_fp, misaligned_reads):
  csv_fp.writelines((', '.join(map(str, line)) + '\n' for line in misaligned_reads))


def main(bam_fp, csv_fp, json_fp, window, block_size):
  cnt = 0
  gen = analyze_bam(bam_fp, window=window, block_size=block_size)
  write_csv_header(csv_fp)
  for mis_reads, tot_reads_ctr, bad_reads_ctr in gen:
    write_csv(csv_fp, mis_reads)
    cnt += len(mis_reads)
    logger.debug('{:d} mis-aligned reads'.format(cnt))
  json.dump({"read_counts": {str(k): v for k,v in tot_reads_ctr.iteritems()},
             "bad_read_counts": {str(k): v for k,v in bad_reads_ctr.iteritems()}},
            json_fp, indent=2)


def cli():
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  with pysam.Samfile(args['--inbam']) as bam_in_fp, \
       open(args['--out_prefix'] + '.csv', 'w') as csv_out_fp, \
       open(args['--out_prefix'] + '.json', 'w') as json_out_fp:
    main(bam_fp=bam_in_fp, csv_fp=csv_out_fp, json_fp=json_out_fp, window=int(args['--window']), block_size=int(args['--block_size']))

if __name__ == "__main__":
  cli()
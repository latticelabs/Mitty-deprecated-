"""Given a bam file containing simulated reads aligned by a tool produce a table of mis-aligned reads. The table has the
following columns

seq,  correct_chrom, correct_pos, aligned_chrom, aligned_pos, mapping_qual, pair mapped

Commandline::

  Usage:
    checkbam  --inbam=INBAM  --fout=FOUT  [--block_size=BS] [-v]

  Options:
    --inbam=INBAM           Input bam file name of reads
    --fout=FOUT             Name of output file
    --block_size=BS         Number of reads to process at a time [default: 1000000]
    -v                      Dump detailed logger messages
"""
__version__ = '1.0.0'

import pysam
import docopt
from mitty.vcf2reads import interpret_read_qname
import logging
logger = logging.getLogger(__name__)


def analyze_bam(bam_fp, block_size=100000):
  """Given a read iterator from a BAM file, process each read and return a list of lists suitable for writing to csv
  file."""
  keep_reading = True
  while keep_reading:
    misaligned_reads = [None] * block_size
    current_count = block_size
    while current_count > 0:
      try:
        read = bam_fp.next()
        correct_chrom_no, _, correct_pos, correct_cigar = interpret_read_qname(read.qname, read.is_read2)
        if correct_chrom_no == read.tid + 1 and correct_pos == read.pos + 1:  #BAM is zero indexed
          continue
        misaligned_reads[block_size - current_count] = [read.seq, correct_chrom_no, correct_pos, read.tid + 1, read.pos, read.mapq, False]
        current_count -= 1
      except StopIteration:
        keep_reading = False
        break
    if current_count > 0:
      del misaligned_reads[block_size - current_count:]

    yield misaligned_reads


def write_csv_header(csv_fp):
  csv_fp.write('seq,  correct_chrom, correct_pos, aligned_chrom, aligned_pos, mapping_qual, pair mapped\n')


def write_csv(csv_fp, misaligned_reads):
  for line in misaligned_reads:
    csv_fp.write(','.join(map(str, line)) + '\n')


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  with pysam.Samfile(args['--inbam']) as bam_in_fp, open(args['--fout'], 'w') as csv_out_fp:
    cnt = 0
    gen = analyze_bam(bam_in_fp, block_size=int(args['--block_size']))
    write_csv_header(csv_out_fp)
    for mis_reads in gen:
      write_csv(csv_out_fp, mis_reads)
      cnt += len(mis_reads)
      logger.debug('{:d} mis-aligned reads'.format(cnt))
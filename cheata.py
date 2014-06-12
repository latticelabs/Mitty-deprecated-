"""This script reads in a BAM file created by `reads.py` and aligns the reads by using the coordinates stored in the
seq id string. This aligned file can be read in using a visualizer to debug the simulation chain. Using the `split`
option treats the input file as an BAM produced by an aligner and splits the reads into two files: correctly aligned
and incorrectly aligned (_correct.bam, _wrong.bam).

Usage:
cheata  --inbam=INBAM  --outbam=OUTBAM  --heada=HD  [-v]
cheata split --inbam=INBAM  [-v]

Options:
  --inbam=INBAM           Input (unaligned) bam file name of reads from reads.py
  --outbam=OUTBAM         Output (perfectly aligned) bam file name
  --heada=HD              Heada file produced by converta. Contains sequence name and len. These are put into the
                          cheat aligned BAM file header as they are needed by some tools such as Tablet
  -v                      Dump detailed logger messages

Notes:
1. Recall that the seq id of each read is the string 'rN:POS1:CIGAR1:POS2:CIGAR2'. Unpaired reads have no POS2, CIGAR2
2. As of v 0.2.0 there is a quirk with pysam writing BAM files where alignments have some hidden problem when written as
BAM from cheata. My work around is to save as SAM and then convert to BAM which seems to not create any problems.
"""
__version__ = '0.2.0'

import os  # Needed for filesize (for sequence length) and messing with filenames
import pysam  # Needed to read/write BAM files
import docopt
import logging
logger = logging.getLogger(__name__)


def sort_and_index_bam(bamfile):
  """Do the filename gymnastics required to end up with a sorted, indexed, bam file."""
  # samtools sort adds a '.bam' to the end of the file name.
  pysam.sort(bamfile, os.path.splitext(bamfile)[0])
  pysam.index(bamfile)


def align(inbam, outbam, seq_name, seq_len):
  in_bamfile = pysam.Samfile(inbam, 'rb')
  out_hdr = in_bamfile.header
  out_hdr['SQ'][0]['SN'] = seq_name.split(' ')[0].strip()  # Some programs, like tablet, can't handle seq names with spaces
  out_hdr['SQ'][0]['LN'] = int(seq_len)
  out_bamfile = pysam.Samfile(outbam, 'wb', header=out_hdr)

  cnt = 0
  blk = 0
  for read in in_bamfile:
    parts = read.qname.split(':')
    read.mapq = 100  # It's better to set this
    read.pos = int(parts[1])-1  # 0-indexed
    read.cigarstring = parts[2]
    if read.flag & 0x01:  # Paired reads
      if read.flag & 0x80:  # Second end (mate)
        read.pos = int(parts[3])-1
        read.cigarstring = parts[4]
    out_bamfile.write(read)

    blk += 1
    if blk == 100000:
      cnt += blk
      logger.debug('{:d} reads done'.format(cnt))
      blk = 0

  logger.debug('{:d} reads done'.format(cnt))
  out_bamfile.close()
  sort_and_index_bam(outbam)


def split_good_bad(inbam):
  """Split the aligned BAM file into two sets of reads - properly aligned and misaligned - and save them in two bam
  files"""

  correct_bam = os.path.splitext(inbam)[0] + '_correct.bam'
  wrong_bam = os.path.splitext(inbam)[0] + '_wrong.bam'

  in_bamfile = pysam.Samfile(inbam, 'rb')

  out_hdr = in_bamfile.header
  correct_bamfile = pysam.Samfile(correct_bam, 'wb', header=out_hdr)
  wrong_bamfile = pysam.Samfile(wrong_bam, 'wb', header=out_hdr)

  cnt = 0
  blk = 0
  good_cnt = 0
  bad_cnt = 0
  for read in in_bamfile:
    cheat_answer = read.qname.split(':')
    this_correct_pos = int(cheat_answer[1])
    if read.flag & 0x01:  # Paired reads
      if read.flag & 0x80:  # Second end (mate)
        this_correct_pos = int(cheat_answer[3])

    computed_pos = read.pos + 1  # PySam uses 0 indexing ...

    if this_correct_pos != computed_pos:  # We go into the bad pile
      wrong_bamfile.write(read)
      bad_cnt += 1
    else:
      correct_bamfile.write(read)
      good_cnt += 1
    blk += 1
    if blk == 100000:
      cnt += blk
      logger.debug('{:d} reads done'.format(cnt))
      blk = 0

  logger.debug('{:d} reads done'.format(cnt))
  logger.debug('Found {:d} good reads, {:d} bad reads'.format(good_cnt, bad_cnt))

  correct_bamfile.close()
  wrong_bamfile.close()
  sort_and_index_bam(correct_bam)
  sort_and_index_bam(wrong_bam)


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  if args['split']:  # We are in split_good_bad mode
    split_good_bad(args['--inbam'])
  else:  # We are in align mode
    seq_name, seq_len = open(args['--heada']).readlines()
    align(args['--inbam'], args['--outbam'], seq_name, seq_len)


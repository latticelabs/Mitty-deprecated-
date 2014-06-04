"""This script reads in a BAM file created by `reads.py` and aligns the reads by using the coordinates stored in the
seq id string. This aligned file can be read in using a visualizer to debug the simulation chain. Using the `check`
option treats the input file as an BAM produced by an aligner and the check file checks the alignment and produces
some statistics of alignment performace in the check file.

Usage:
cheata  --inbam=INBAM  --outbam=OUTBAM  --heada=HD  [-v]
cheata check --inbam=INBAM  --checkfile=CHECK  [-v]

Options:
  --inbam=INBAM           Input (unaligned) bam file name of reads from reads.py
  --outbam=OUTBAM         Output (perfectly aligned) bam file name
  --heada=HD              Heada file produced by converta. Contains sequence name and len. These are put into the
                          cheat aligned BAM file header as they are needed by some tools such as Tablet
  --checkfile=CHECK       Output file containing alignment analysis statistics
  -v                      Dump detailed logger messages

Notes:
1. Recall that the seq id of each read is the string 'rN:POS1:CIGAR1:POS2:CIGAR2'. Unpaired reads have no POS2, CIGAR2
2. As of v 0.2.0 there is a quirk with pysam writing BAM files where alignments have some hidden problem when written as
BAM from cheata. My work around is to save as SAM and then convert to BAM which seems to not create any problems.
"""
__version__ = '0.2.0'

import tempfile  # Needed because we sort the original alignment and then index it
import os  # Needed for removing the extra .bam samtools sort adds to the name and for filesize (for sequence length)
import pysam  # Needed to read/write BAM files
import cPickle
import docopt
import logging
logger = logging.getLogger(__name__)


def align(inbam, outbam, seq_name, seq_len):
  in_bamfile = pysam.Samfile(inbam, 'rb')

  # We save the reads first to a temporary file, sort them and save the sorted aligned reads under the name we want
  tf_h, tf_name = tempfile.mkstemp()
  os.close(tf_h)

  out_hdr = in_bamfile.header
  out_hdr['SQ'][0]['SN'] = seq_name.split(' ')[0].strip()  # Some programs, like tablet, can't handle seq names with spaces
  out_hdr['SQ'][0]['LN'] = int(seq_len)
  out_bamfile = pysam.Samfile(tf_name, 'wb', header=out_hdr)

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

  pysam.sort(tf_name, outbam)
  os.rename(outbam + '.bam', outbam)
  # samtools sort adds a '.bam' to the end of the file name. This rename reverses that. We chose to rename because we
  # like not to make any assumptions about what extension the user wants for their file. For example, we do not assume
  # that the user's file will end in .bam and simply try and strip that.
  pysam.index(outbam)
  os.remove(tf_name)  # Clean up after ourselves


def check(inbam, checkfile):
  in_bamfile = pysam.Samfile(inbam, 'rb')

  read_id = []
  correct_pos = []
  aligned_pos = []
  bad_alignment_map_score = []  # These scores match with the previous arrays

  good_alignment_map_score = []

  cnt = 0
  blk = 0

  for read in in_bamfile:
    cheat_answer = read.qname.split(':')
    this_correct_pos = int(cheat_answer[1])
    if read.flag & 0x01:  # Paired reads
      if read.flag & 0x80:  # Second end (mate)
        this_correct_pos = int(cheat_answer[3])

    computed_pos = read.pos + 1  # PySam uses 0 indexing ...

    if this_correct_pos != computed_pos:
      read_id.append(read.qname)
      correct_pos.append(this_correct_pos)
      aligned_pos.append(computed_pos)
      bad_alignment_map_score.append(read.mapq)
    else:
      good_alignment_map_score.append(read.mapq)

    blk += 1
    if blk == 100000:
      cnt += blk
      logger.debug('{:d} reads done'.format(cnt))
      blk = 0

  logger.debug('{:d} reads done'.format(cnt))

  cPickle.dump({
    'comment': "correct_pos is the correct position of the misaligned read, aligned_pos is the pos the aligner placed the read at.",
    'seq': in_bamfile.header['SQ'][0]['SN'],  # Name of the reference sequence
    'len': int(in_bamfile.header['SQ'][0]['LN']),  # Length of the reference sequence
    'read_ids': read_id,
    'correct_pos': correct_pos,
    'aligned_pos': aligned_pos,
    'bad_alignment_map_score': bad_alignment_map_score,
    'good_alignment_map_score': good_alignment_map_score
  }, open(checkfile, 'w'), protocol=cPickle.HIGHEST_PROTOCOL)



if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  if args['check']:  # We are in check mode
    check(args['--inbam'], args['--checkfile'])
  else:  # We are in align mode
    seq_name, seq_len = open(args['--heada']).readlines()
    align(args['--inbam'], args['--outbam'], seq_name, seq_len)


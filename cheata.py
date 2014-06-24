"""In default mode, this script reads in a BAM file created by `reads.py` and aligns the reads by using the coordinates
stored in the seq id string. This aligned file can be read in using a visualizer to debug the simulation chain.

Using the `split` option treats the input file as an BAM produced by an aligner and splits the reads into
correctly aligned, incorrectly aligned and unmapped (_correct.bam, _wrong.bam, _unmapped.bam).

Usage:
cheata  --inbam=INBAM  --outbam=OUTBAM  --heada=HD  [-v]
cheata split --inbam=INBAM  --wg=WG [-v]

Options:
  --inbam=INBAM           Input (unaligned) bam file name of reads from reads.py
  --outbam=OUTBAM         Output (perfectly aligned) bam file name
  --wg=WG                 Whole genome file. Needed for the index, to match BAM seq ids to chrom_key in read qname
  --heada=HD              Heada file produced by converta. Contains sequence name and len. These are put into the
                          cheat aligned BAM file header as they are needed by some tools such as Tablet
  -v                      Dump detailed logger messages

Notes:
1. Recall that the seq id of each read is the string 'chrom:copy:rN:POS1:CIGAR1:POS2:CIGAR2'.
   Unpaired reads have no POS2, CIGAR2
2. As of v 0.2.0 there is a quirk with pysam writing BAM files where alignments have some hidden problem when written as
   BAM from cheata. My work around is to save as SAM and then convert to BAM which seems to not create any problems.
3. The `split` mode enables us to break the larger aligned .bam file into smaller files. Since we are often interested
   in just the errors, we can look at just the _wrong.bam which is often much smaller.
"""
__version__ = '0.2.0'

import os  # Needed for filesize (for sequence length) and messing with filenames
import pysam  # Needed to read/write BAM files
import docopt
import genome
import reads
import logging
logger = logging.getLogger(__name__)


def sort_and_index_bam(bamfile):
  """Do the filename gymnastics required to end up with a sorted, indexed, bam file."""
  # samtools sort adds a '.bam' to the end of the file name.
  pysam.sort(bamfile, os.path.splitext(bamfile)[0])
  pysam.index(bamfile)


def align(inbam, outbam, seq_name, seq_len):
  logger.error('Not updated yet')
  raise RuntimeError

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


def split(wg, in_bam, correct_bam, wrong_bam, unmapped_bam):
  """Split the aligned BAM file into three sets of reads - properly aligned, misaligned and unmapped - and save them in
  separate bam files
  """
  #The rname in the aligned read corresponds to the sequence in the list in_bam.header['SQ']
  #We need to match that up with the chromosome code we have in the wg file
  seq_name_list = [sq['SN'] for sq in in_bam.header['SQ']]
  rev_idx = wg.reverse_index
  cnt = 0
  blk = 0
  good_cnt, bad_cnt, unmapped_cnt = 0, 0, 0
  for read in in_bam:
    if read.flag & 0x4:  # Unmapped
      unmapped_bam.write(read)
      unmapped_cnt += 1
    else:
      correct_chrom, correct_pos, correct_cigar = reads.interpret_read_qname(read)
      aligned_chrom = rev_idx[seq_name_list[read.rname]]
      aligned_pos = read.pos + 1  # PySam uses 0 indexing ...
      if correct_chrom == aligned_chrom and correct_pos == aligned_pos:
        correct_bam.write(read)
        good_cnt += 1
      else:
        wrong_bam.write(read)
        bad_cnt += 1
    blk += 1
    if blk == 100000:
      cnt += blk
      logger.debug('{:d} reads done'.format(cnt))
      blk = 0

  logger.debug('{:d} reads done'.format(cnt))
  logger.debug('Found {:d} good reads, {:d} bad reads, {:d} unmapped reads'.format(good_cnt, bad_cnt, unmapped_cnt))


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  if args['split']:  # We are in split_good_bad mode
    inbam_name = args['--inbam']
    correct_bam_name = os.path.splitext(inbam_name)[0] + '_correct.bam'
    wrong_bam_name = os.path.splitext(inbam_name)[0] + '_wrong.bam'
    unmapped_bam_name = os.path.splitext(inbam_name)[0] + '_unmapped.bam'
    with genome.WholeGenome(fname=args['--wg']) as wg:
      in_bam_f=pysam.Samfile(inbam_name, 'rb')
      split(wg, in_bam=in_bam_f,
            correct_bam=pysam.Samfile(correct_bam_name, 'wb', header=in_bam_f.header),
            wrong_bam=pysam.Samfile(wrong_bam_name, 'wb', header=in_bam_f.header),
            unmapped_bam=pysam.Samfile(unmapped_bam_name, 'wb', header=in_bam_f.header))
  else:  # We are in align mode
    seq_name, seq_len = open(args['--heada']).readlines()
    align(args['--inbam'], args['--outbam'], seq_name, seq_len)


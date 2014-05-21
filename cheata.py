"""This script reads in a BAM file created by `reads.py` and 'aligns' the reads by using the coordinates stored in the
seq id string. This 'aligned' file can be read in using a visualizer to debug the simulation chain.

Usage:
cheata --inbam=BAM  --outbam=BAM  --refname=REFNAME  --reflen=REFLEN  [-v]

Options:
  --inbam=INBAM           Input (unaligned) bam file name of reads from reads.py
  --outbam=OUTBAM         Output (perfectly aligned) bam file name
  --refname=REFNAME       Reference sequence name
  --reflen=REFLEN         Reference sequence lengths
                          (These are needed by other tools that look at the BAM file)
  -v                      Dump detailed logger messages

Notes:
1. Recall that the seq id of each read is the string 'rN:S1:S2' where N is the number of the read,
   S1 the start of the first read and S2 the start of the mate pair. Unpaired reads have no S2
2. As of v 0.2.0 there is a quirk with pysam writing BAM files where alignments have some hidded problem when written as
BAM from cheatah. My wor around is to save as SAM and then convert to BAM which seems to not create any problems.
"""
__version__ = '0.2.0'

import tempfile  # Needed because we sort the original alignment and then index it
import os  # Needed for removing the extra .bam samtools sort adds to the name and for filesize (for sequence length)
import pysam  # Needed to read/write BAM files
import docopt


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  in_bamfile = pysam.Samfile(args['--inbam'], 'rb')

  # seq_name = open(args['--refseq'] + '.heada', 'r').readline()
  # seq_len = os.path.getsize(args['--refseq'])

  seq_name = args['--refname']
  seq_len = int(args['--reflen'])

  out_hdr = in_bamfile.header
  out_hdr['SQ'][0]['SN'] = seq_name.split(' ')[0]  # Some programs, like tablet, can't handle seq names with spaces
  out_hdr['SQ'][0]['LN'] = seq_len

  # We save the reads first to a temporary file, sort them and save the sorted aligned reads under the name we want
  tf_h, tf_name = tempfile.mkstemp()
  os.close(tf_h)
  out_bamfile = pysam.Samfile(tf_name, 'wb', header=out_hdr)

  cnt = 0
  blk = 0
  for read in in_bamfile:
    parts = read.qname.split(':')
    aligned_read = pysam.AlignedRead()
    aligned_read.qname = read.qname
    aligned_read.seq = read.seq
    aligned_read.qual = read.qual
    aligned_read.mapq = 100  # It's better to set this
    aligned_read.pos = int(parts[1])-1  # 0-indexed
    aligned_read.cigarstring = parts[2]
    aligned_read.flag = read.flag
    if read.flag & 0x01:  # Paired reads
      if read.flag & 0x80:  # Second end (mate)
        aligned_read.pos = int(parts[3])-1
        aligned_read.cigarstring = parts[4]
    out_bamfile.write(aligned_read)

    cnt += 1
    blk += 1
    if blk == 100000:
      print '{:d} reads done'.format(cnt)
      blk = 0

  print '{:d} reads done'.format(cnt)
  out_bamfile.close()

  pysam.sort(tf_name, args['--outbam'])
  os.rename(args['--outbam'] + '.bam', args['--outbam'])
  # samtools sort adds a '.bam' to the end of the file name. This rename reverses that. We chose to rename because we
  # like not to make any assumptions about what extension the user wants for their file. For example, we do not assume
  # that the user's file will end in .bam and simply try and strip that.
  pysam.index(args['--outbam'])
  os.remove(tf_name)  # Clean up after ourselves


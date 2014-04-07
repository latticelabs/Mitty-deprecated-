"""This script reads in a BAM file created by `reads.py` and 'aligns' the reads by using the coordinates stored in the
seq id string. This 'aligned' file can be read in using a visualizer to debug the simulation chain.

Usage:
cheata --inbam=BAM  --outbam=BAM  [-v]

Options:
  --inbam=INBAM           Input (unaligned) bam file name of reads from reads.py
  --outbam=OUTBAM         Output (perfectly aligned) bam file name
  -v                      Dump detailed logger messages

Notes:
1. Recall that the seq id of each read is the string 'rN:S1:S2' where N is the number of the read,
   S1 the start of the first read and S2 the start of the mate pair. Unpaired reads have no S2
"""
__version__ = '0.1.0'

import pysam  # Needed to read/write BAM files
import docopt


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  in_bamfile = pysam.Samfile(args['--inbam'], 'rb')
  out_bamfile = pysam.Samfile(args['--outbam'], 'wb', header=in_bamfile.header)

  iter = in_bamfile.fetch()
  for x in iter:
    parts = x.qname.split(':')
    x.flag |= 0x2  # We've mapped this segment
    if x.flag & 0x01:
      if x.flag & 0x40:
        x.pos = int(parts[1])-1
      else:
        x.pos = int(parts[2])-1
    else:
      x.pos = int(parts[1])-1
    out_bamfile.write(x)
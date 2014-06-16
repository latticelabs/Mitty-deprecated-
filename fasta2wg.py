"""This converts whole genome fasta files e.g. hg38.fa.gz - "Soft-masked" assembly sequence in one file
released by the GRC into a .wg.gz file.

Usage:
fasta2wg  --fa=FA  --wg=WG  --chrom_count=CC  [--species=SP] [-v]

Options:
  --fa=FA           Input .fa.gz file
  --wg=WG           Name of .wg.gz file to save to (Saved as gzipped .wg file)
  --chrom_count=CC  How many chromosomes in the file
  --species=SP      What string would we like at appear in the whole genome file [default: Test]
  -v                Be verbose when you do things

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
import gzip
import docopt
import genome
import logging
logger = logging.getLogger(__name__)

__version__ = '0.2.0'


def fasta_read_generator(fasta_fname):
  """Given a fasta filename return a generator that yields us sequences and sequence ids from the file"""
  with gzip.open(fasta_fname, 'rb') as fasta_fp:
    read_a_byte = fasta_fp.read
    reading_a_seq = False
    seq_id = None
    seq = []
    while True:
      byte = read_a_byte(1)
      if not byte or byte == '>': # Either file is done, or a new seq starting
        if reading_a_seq:
          fasta_fp.seek(-1, 1)  # Pretend this never happened next time round
          reading_a_seq = False
          yield seq_id, ''.join(seq)
        if byte == '>':  # New sequence
          reading_a_seq = True
          seq_id = fasta_fp.readline()[:-1]
          seq = []
        if not byte:  # Done with file
          break
      elif byte != '\n':  # Part of the sequence
        seq += byte.upper()

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  with genome.WholeGenome(fname=args['--wg'], species=args['--species'],
                          chrom_count=int(args['--chrom_count'])) as wg:
    for n, (this_seq_id, this_seq) in enumerate(fasta_read_generator(args['--fa'])):
      if wg.insert_seq(this_seq, chrom_no=n + 1, chrom_cpy=1, seq_id=this_seq_id):
        print 'Inserted chromosome {:d}'.format(n + 1)

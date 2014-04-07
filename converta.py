"""Script to convert from fasta to smalla files

Usage:
converta <fasta> <smalla> [--block_size=BS]

Options:
  fasta              Input file name
  smalla             Output file name
  --block_size=BS    Block size [default: 1000000]

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
import docopt
import seqio

__version__ = '0.1.0'


def fasta_to_smalla(fasta_fname, smalla_fname, block_size=1000):
  with open(smalla_fname, 'w') as f:
    for s in seqio.block_read_fasta(fname=fasta_fname, block_size=block_size):
      f.write(s)


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)
  fasta_to_smalla(args['<fasta>'], args['<smalla>'], int(args['--block_size']))

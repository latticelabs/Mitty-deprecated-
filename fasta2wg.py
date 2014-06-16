"""This converts whole genome fasta files e.g. hg38.fa.gz - "Soft-masked" assembly sequence in one file
released by the GRC into a .wg.gz file.

Usage:
fasta2wg  --index=IDX  --wg=WG  [-v]

Options:
  --index=IDX       Index file (.json) listing fasta files to be inserted into whole genome
  --wg=WG           Name of .wg.gz file to save to (Saved as gzipped .wg file)
  -v                Be verbose when you do things

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
import json
import gzip
import docopt
import genome
import logging
logger = logging.getLogger(__name__)

__version__ = '0.2.0'


def read_single_seq_fasta(fasta_fname):
  """Given a fasta filename with one sequence read it."""
  with gzip.open(fasta_fname, 'rb') as fasta_fp:
    seq_id = fasta_fp.readline()[1:-1]
    seq = fasta_fp.read().replace('\n','').upper()
  return seq_id, seq

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  idx = json.load(open(args['--index'],'r'))
  with genome.WholeGenome(fname=args['--wg'], species=idx['header']['species'].encode('ascii'),
                          chrom_count=idx['header']['chromosome count']) as wg:
    for k, fasta_fname in idx['chromosomes'].iteritems():

      chrom_no, chrom_cpy = int(k.split(':')[0]), int(k.split(':')[1])
      this_seq_id, this_seq = read_single_seq_fasta(fasta_fname)
      if wg.insert_seq(this_seq, chrom_no=chrom_no, chrom_cpy=chrom_cpy, seq_id=this_seq_id):
        print 'Inserted chromosome {:d}, copy {:d} ({:s})'.format(chrom_no, chrom_cpy, this_seq_id)

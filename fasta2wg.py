"""Given a list of fasta files, each containing one fasta sequence we compact this into a whole genome file. The
list of files is given in the index file (.json)

Usage:
fasta2wg  --index=IDX  --wg=WG  [--fa=FA] [-v]
fasta2wg  explain

Options:
  --index=IDX       Index file (.json) listing fasta files to be inserted into whole genome
  --wg=WG           Name of .wg.gz file to save to (Saved as gzipped .wg file)
  --fa=FA           If set, will also dump a fa.gz file (by simply concatenating the files together) for use by BWA
  -v                Be verbose when you do things
  explain           Print index file format and exit

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
import h5py
import numpy
import json
import gzip
import docopt
import logging
logger = logging.getLogger(__name__)

__version__ = '0.2.0'

__explain__ = """
Index file example:

{
    "header": {
        "species": "Test Chimera",
    },
    "chromosomes": {
        "1:1": "Data/porcine_circovirus.fa.gz",
        "2:1": "Data/adenovirus.fa.gz",
        "1:2": "Data/altered_porcine.fa.gz",
        "2:2": "Data/herpes.fa.gz",
        "3:1": "Data/parvovirus.fa.gz"
    }
}

f['sequence/1/1'][:30].tostring()
"""


def save_genome_to_hdf5(index, h5_fp):
  """Given a fasta filename with one sequence, read it."""
  h5_fp.attrs['species'] = index['header']['species']
  grp = h5_fp.create_group('sequence')
  for k, fasta_fname in index['chromosomes'].iteritems():
    chrom_no, chrom_cpy = int(k.split(':')[0]), int(k.split(':')[1])
    with gzip.open(fasta_fname, 'rb') as fasta_fp:
      seq_id = fasta_fp.readline()[1:-1]
      dset = h5_fp.create_dataset('{:s}/{:d}/{:d}'.format(grp.name, chrom_no, chrom_cpy),
                                  data=numpy.fromstring(fasta_fp.read().replace('\n', '').upper(), dtype='u1'))
      dset.attrs['seq id'] = seq_id
      logger.debug('Inserted chromosome {:d}, copy {:d} ({:s})'.format(chrom_no, chrom_cpy, seq_id))


def concatenate_fasta(file_list, fasta_out):
  """A wrapper around cat and gzip.
  1. You can concatenate gzipped files and they will work!
  2. Concatenating chains the strings together, but we need a newline after each sequence
  """
  import os, subprocess
  # This is our newline after each file
  with gzip.open(fasta_out + '.sp', 'w') as fp:
    fp.write('\n')
  f_list = ' {:s} '.format(fasta_out + '.sp').join(file_list)
  #TODO watch out for overly long command lines here
  _ = subprocess.call('cat {:s} > {:s}'.format(f_list, fasta_out), shell=True)  # We live dangerously
  os.remove(fasta_out + '.sp')  # Clean up after ourselves

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  if args['explain']:
    print __explain__
    exit(0)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  idx = json.load(open(args['--index'],'r'))
  f = h5py.File(args['--wg'], "w")
  save_genome_to_hdf5(idx, f)
  f.close()

  if args['--fa']:
    concatenate_fasta(idx['chromosomes'].values(), args['--fa'])
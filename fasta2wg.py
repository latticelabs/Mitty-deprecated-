"""Given a list of fasta files, each containing one fasta sequence we compact this into a whole genome file.
list of files is given in the index file (.json)

Usage:
fasta2wg  --index=IDX  --wg=WG  [--fa=FA] [-v]
fasta2wg  describe  --wg=WG
fasta2wg  explain

Options:
  --index=IDX       Index file (.json) listing fasta files to be inserted into whole genome
  --wg=WG           Name of whole genome file
  --fa=FA           If set, will also dump a fa.gz file (by simply concatenating the files together) for use by BWA
  -v                Be verbose when you do things
  describe          Take the indicated whole genome file as input and print a summary of what it contains
  explain           Print index file format and exit

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
import h5py
import os
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
    "chromosomes": [
        ["Data/porcine_circovirus.fa.gz"],
        ["Data/adenovirus.fa.gz"],
        ["Data/altered_porcine.fa.gz"],
        ["Data/herpes.fa.gz"],
        ["Data/parvovirus.fa.gz"]
    ]
}

f['sequence/1/1'][:30].tostring()
"""


def save_genome_to_hdf5(index, h5_fp):
  """Store all the fasta sequences in the appropriate places in the hdf5 file.
  Chromosomes are stored under 'sequence/1' 'sequence/2' ....
  The sequence id (say from NIH) is stored as a group attribute here
  The seq id is taken as the seq id of the first
  Copies of chromosomes are stored as 1,2 ...

  """
  h5_fp.attrs['species'] = index['header']['species'].encode('ascii')  # Ubuntu HDF5 version, issue with unicode
  grp = h5_fp.create_group('sequence')
  for n, fasta_fnames in enumerate(index['chromosomes']):
    chrom_grp = grp.create_group(str(n + 1))
    for m, fasta_fname in enumerate(fasta_fnames):
      with gzip.open(fasta_fname, 'rb') as fasta_fp:
        seq_id = fasta_fp.readline()[1:-1]
        if 'seq_id' not in chrom_grp.attrs:
          chrom_grp.attrs['seq_id'] = seq_id
        h5_fp.create_dataset('{:s}/{:d}/{:d}'.format(grp.name, n + 1, m + 1),
                                    data=numpy.fromstring(fasta_fp.read().replace('\n', '').upper(), dtype='u1'))
        logger.debug('Inserted chromosome {:d}, copy {:d} ({:s})'.format(n + 1, m + 1, seq_id))


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


def describe(h5_fp):
  print 'Species: {:s}'.format(h5_fp.attrs['species'])
  print 'Chromosomes: {:d}'.format(len(h5_fp['sequence']))
  for chrom, v in h5_fp['sequence'].iteritems():
    print 'Chrom {:s} ({:s})'.format(chrom, v.attrs['seq_id'])
    for cpy, seq in v.iteritems():
      mutated = '(mutated)' if 'pos/{:s}/{:s}'.format(chrom, cpy) in h5_fp else '(reference)'
      print '\tCopy {:s} {:d}bp {:s}'.format(cpy, seq.size, mutated)


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

  if args['describe']:
    with h5py.File(args['--wg'], 'r') as fp:
      describe(fp)
  else:
    if not os.path.isdir(os.path.dirname(args['--wg'])):  # Not handling race conditions
        os.makedirs(os.path.dirname(args['--wg']))
    with h5py.File(args['--wg'], 'w') as fp:
      idx = json.load(open(args['--index'], 'r'))
      save_genome_to_hdf5(idx, fp)
    if args['--fa']:
      concatenate_fasta([fname for fname_groups in idx['chromosomes'] for fname in fname_groups], args['--fa'])
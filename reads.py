"""This module contains functions that generate simulated reads. The module can be called
as a script as well. This is useful for creating test data for MGR algorithms/data formats

Usage:
reads [paired] [options]

Options:
  --ref=REF                  Reference chromosome [default: ../../../Data/GRCh38/chr24.fa]
  --fastq=FASTQ              FASTQ output file name (no extension) [default: simulated_reads]
  --seed=SEED                Seed for random number generator [default: 0]
  --read_len=READLEN         Length of reads [default: 100]
  --paired_len=PAIRLEN       Length of whole segment for paired reads [default: 1000]
  --read_count=READCNT       Number of reads [default: 1000]
  --comment=COMMENT          User comment [default: No comment]

Notes:
1. The output is a FASTQ file (or two FASTQ files for paired reads).
2. The quality scores are in Sanger/Phred format
3. The seq id contains a special prefix for this tool, a sequential read number and coordinates for the read
4. Any annotations the user wishes to make (plus some tool info) are stored in a sidecar file with the same name as
   the FASTQ file with .info added to the end
5. For paired reads two FASTQ files are created, numbered _1.fastq and _2.fastq
6. From informal testing the major time bottle neck is in writing out the fastq file(s)
"""
#__version__ = '0.1.0' #Fixed read lengths. Uniform coverage
__version__ = '0.2.0'  #Paired end. Fixed read lengths. Uniform coverage

import docopt
import numpy
from Bio import SeqIO  # Need for loading reference sequence
import logging

logger = logging.getLogger(__name__)


def generate_reads(reference, read_len=100, paired_len=1000, num_reads=1000, paired=False, seed=0):
  """Given a reference sequence generate reads of a given length.

  Inputs
    reference  - string(like) containing DNA sequence
    read_len   - the length of the reads
    num_reads  - the number of reads we will generate
    seed       - the seed for the random number generator

  Outputs
    A list of one or two dictionaries with keys
      'sequences'  - a list of sequence strings (or any other indexable dtype) that represent base letters
      'indexes'    - a list of tuples indicating the start and stop of the reads
      'qualities'  - a list of lists of integers between 0 and 96 representing base call quality (phred score)
    If we have single reads there is only one dictionary in the list, if we have paired reads then there are two
    dictionaries in the list
  """
  rng = numpy.random.RandomState(seed)  # Initialize the numpy RNG
  logger.debug('Starting to generate reads')
  if paired:
    read_starts = rng.randint(0, len(reference) - paired_len, num_reads)
    read1 = {
      'sequences': [reference[read_starts[n]:read_starts[n] + read_len] for n in range(num_reads)],
      'indexes': [(st, st + read_len) for st in read_starts],
      'qualities': [[96] * read_len] * num_reads
    }
    read2 = {
      'sequences': [reference[read_starts[n] + paired_len - read_len:read_starts[n] + paired_len] for n in
                    range(num_reads)],
      'indexes': [(st + paired_len - read_len, st + paired_len) for st in read_starts],
      'qualities': [[96] * read_len] * num_reads
    }
    reads = [read1, read2]
  else:
    read_starts = rng.randint(0, len(reference) - read_len, num_reads)
    read1 = {
      'sequences': [reference[read_starts[n]:read_starts[n] + read_len] for n in range(num_reads)],
      'indexes': [(st, st + read_len) for st in read_starts],
      'qualities': [[96] * read_len] * num_reads
    }
    reads = [read1]
  logger.debug('Finished generating reads')
  return reads


def write_fastq(reads, file_handles, seq_id_prefix='SBG_sim'):
  """
  Given a list of sequences and their quality write the read data to an already opened text file.

  Inputs:
    reads        -  A list of one or two dictionaries with keys
                  'sequences'  - a list of sequence strings (or any other indexable dtype) that represent base letters
                  'indexes'    - a list of tuples indicating the start and stop of the reads
                  'qualities'  - a list of lists of integers between 0 and 96 representing base call phred score
                  If we have single reads there is only one dictionary in the list, if we have paired reads then there
                  are two dictionaries in the list
    file_handles - a list of file handles. The output will be appended to this file.
                   If we have paired reads there will be two handles, otherwise there will be one
    seq_id_prefix - a string that we attach to the read id

  Notes:
  1. The seq id string contains the coordinates of the bases so we can debug alignment programs
  2. The quality is written out using the sanger format (Phred+33)
  """
  for these_reads, file_handle in zip(reads, file_handles):
    logger.debug('Starting to write fastq')
    for n, (sequence, idx, quality) in enumerate(
        zip(these_reads['sequences'], these_reads['indexes'], these_reads['qualities'])):
      file_handle.write('@{:s} {:12d} {:s}\n{:s}\n+\n{:s}\n'.format(seq_id_prefix, n, str(idx),
                                                                    sequence, ''.join([chr(q) for q in quality])))
    logger.debug('Finished writing fastq')


def write_sidecar(args, file_handle):
  """
  Write parameters into a sidecar file
  Inputs:
    args        - the program arguments as parsed by docopts
    file_handle - handle of an opened text file. The output will be appended to this file.
  """
  file_handle.write('SBG Read simulator v{:s}\n'.format(__version__))
  for k,v in args.iteritems():
    file_handle.write('{:s}: {:s}\n'.format(k, str(v)))


def main(args):
  reference = SeqIO.parse(args['--ref'],'fasta').next().seq #We always have a reference
  reads = generate_reads(reference, read_len=int(args['--read_len']), paired_len=int(args['--paired_len']),
                         num_reads=int(args['--read_count']), paired=args['paired'], seed=int(args['--seed']))
  if len(reads) == 2:
    f1 = open(args['--fastq'] + '_1.fastq', 'w')
    f2 = open(args['--fastq'] + '_2.fastq', 'w')
    file_handles = [f1, f2]
  else:
    file_handles = [open(args['--fastq'] + '.fastq', 'w')]
  write_fastq(reads, file_handles, seq_id_prefix='SBG_READ_SIM {:s} '.format(__version__))

  with open(args['--fastq'] + '.info','w') as file_handle:
    write_sidecar(args, file_handle)

if __name__ == "__main__":
  logging.basicConfig(level=logging.DEBUG)
  arguments = docopt.docopt(__doc__, version=__version__)
  main(arguments)
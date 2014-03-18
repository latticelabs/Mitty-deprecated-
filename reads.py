"""This module contains functions that generate simulated reads. The module can be called
as a script as well. This is useful for creating test data for MGR algorithms/data formats

Usage:
reads [--ref <REF>] [--fastq <FQ>] [--read_len <LEN>] [--read_count <CNT>] [--comment <CMT>]

Options:
  --ref=<REF>                Reference chromosome [default: ../../../Data/GRCh38/chr24.fa]
  --fastq=<FQ>               FASTQ output file name [default: simulated_reads.fastq]
  --read_len=<LEN>           Length of reads [default: 10]
  --read_count=<CNT>         Number of reads [default: 10]
  --comment=<CMT>            User comment [default: No comment]

If the variant call file name is not supplied we generate the NULL model: reads from the reference genome.

Notes:
1. The output is a FASTQ file.
2. The quality scores are in Sanger/Phred format
3. The seq id contains a special prefix for this tool and a sequential read number
4. Any annotations the user wishes to make (plus some tool info) are stored in a sidecar file with the same prefix as
   the FASTQ file ending in .info

"""
__version__ = '0.1.0' #Null model. Fixed read lengths. Uniform coverage

import docopt, numpy, os
from Bio import SeqIO #Need for loading reference sequence

def read_vcf():
  """For future expansion. Currently we will only generate the NULL model"""


def generate_null_reads(reference, read_len=500, num_reads=1000):
  """Given a reference sequence generate reads of a given length.

  Inputs
    reference  - string(like) containing DNA sequence
    read_len   - the length of the reads
    num_reads  - the number of reads we will generate

  Outputs
    sequences  - a list of sequence strings (or any other indexable dtype) that represent base letters
    indexes    - a list of tuples indicating the start and stop of the reads
    qualities  - a list of lists of integers between 0 and 96 representing base call quality (phred score)
  """
  read_starts = numpy.random.randint(0, len(reference) - read_len, num_reads)
  sequences = [reference[read_starts[n]:read_starts[n]+read_len] for n in range(num_reads)]
  qualities = [[96]*read_len]*num_reads
  return sequences, [(st, st+read_len) for st in read_starts], qualities

def write_fastq(sequences, indexes, qualities, file_handle, seq_id_prefix='SBG_sim'):
  """
  Given a list of sequences and their quality write the read data to an already opened text file.

  Inputs:
    sequences   - a list of sequence strings (or any other indexable dtype) that represent base letters
    indexes     - a list of tuples indicating the start and stop of the reads
    qualities   - a list of lists of integers between 0 and 96 representing base call quality/phred score
    file_handle - handle of an opened text file. The output will be appended to this file.

  Notes:
  The quality is written out using the sanger format (Phred+33)
  """
  for n,(sequence,idx,quality) in enumerate(zip(sequences, indexes, qualities)):
    file_handle.write('@{:s} {:12d} {:s}\n{:s}\n+\n{:s}\n'.format(seq_id_prefix,n, str(idx),sequence,''.join([chr(q) for q in quality])))

def write_sidecar(comment, file_handle):
  file_handle.write('SBG Read simulator v{:s}\n User comment'.format(__version__, comment))


def main(args):
  reference = SeqIO.parse(args['--ref'],'fasta').next().seq #We always have a reference
  sequences, indexes, qualities = generate_null_reads(reference,
                                             read_len=int(args['--read_len']), num_reads=int(args['--read_count']))
  with open(args['--fastq'],'w') as file_handle:
    write_fastq(sequences, indexes, qualities, file_handle, seq_id_prefix='SBG_READ_SIM {:s} '.format(__version__))

  with open(args['--fastq'] + '.info','w') as file_handle:
    write_sidecar(args['--comment'], file_handle)

if __name__ == "__main__":
  arguments = docopt.docopt(__doc__, version=__version__)
  main(arguments)
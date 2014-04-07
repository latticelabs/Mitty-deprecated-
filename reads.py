"""This module contains functions that generate simulated reads. The module can be called as a script as well. This is
useful for creating test data for MGR algorithms/data formats

Usage:
reads --paramfile=PFILE [--reads_per_block=BL] [-v]

Options:
  --paramfile=PFILE       Name for parameter file
  --reads_per_block=BL    Generate these many reads at a time (Adjust to machine resources). [default: 10000]
  -v                      Dump detailed logger messages

Notes:
1. The quality scores are in Phred scale (as specified in the SAM spec)
2. Any annotations the user wishes to make (plus the command line arguments and all other parameters used to run the
   sim) are stored in a sidecar file with the same name as the bam file with .info added to the end

"""
#__version__ = '0.1.0'  # Fixed read lengths. Uniform coverage
#__version__ = '0.2.0'  # Paired end. Fixed read lengths. Uniform coverage
__version__ = '0.2.1'  # Plugin system implemented

import imp
import mmap
import docopt
import pysam  # Needed to write BAM files
import logging

logger = logging.getLogger(__name__)


def save_reads_to_bam(bam_file, reads):
  """
  Inputs:
                                 _________ ( seq_str, quality_str, coordinate)
    reads     -  [              /
                  [( ... ), ( ...)],
                  [( ... ), ( ...)], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads
  """
  if len(reads) == 0: return
  paired = True if len(reads[0]) == 2 else False
  for n in range(len(reads)):
    ar = pysam.AlignedRead()
    ar.qname = "r{:d}".format(n)  # Figure out how to store coordinate
    ar.seq = reads[n][0][0]
    ar.qual = reads[n][0][1]
    if paired:
      ar.flag = 0x41  # end1 0x01 flag has to be set to indicate multiple segments
    bam_file.write(ar)

    if paired:
      ar.seq = reads[n][1][0]
      ar.qual = reads[n][1][1]
      ar.flag = 0x81  # end2 0x01 flag has to be set to indicate multiple segments
      bam_file.write(ar)


def write_sidecar(args, params_file, file_handle):
  """
  Write parameters into a sidecar file
  Inputs:
    args        - the program arguments as parsed by docopts
    file_handle - handle of an opened text file. The output will be appended to this file.
  """
  file_handle.write('SBG Read simulator v{:s}\n'.format(__version__))
  file_handle.write('Commandline:\n')
  for k,v in args.iteritems():
    file_handle.write('{:s}: {:s}\n'.format(k, str(v)))
  file_handle.write('Parameters:\n')
  with open(params_file, 'r') as f:
    for line in f.readlines():
      file_handle.write(line)


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__,['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)
  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=logging.DEBUG)

  params = imp.load_source('params', args['--paramfile'], open(args['--paramfile'], 'r'))

  #Load the ref-seq smalla file
  f = open(params.seq_file, 'r+b')
  ref_seq = mmap.mmap(f.fileno(), 0)
  ref_seq_len = len(ref_seq)

  model_fname = 'Plugins/Reads/' + params.model_name + '_plugin.py'  # TODO: join paths properly
  read_model = imp.load_source('readmodel', model_fname, open(model_fname, 'r'))

  if params.stop == -1: params.stop = ref_seq_len
  subsequence_len = params.stop - params.start
  total_reads = int(params.coverage * subsequence_len / float(params.average_read_len))
  reads_per_block = int(args['--reads_per_block'])

  # TODO: drop fastq file as needed
  bam_hdr = {'HD': {'VN': '1.4'},
             'SQ': [{'LN': ref_seq_len, 'SN': params.seq_header,
                     'SP': 'Simulated human, perfect reads, by reads.py {:s}'.format(__version__)}]}
  perfect_bam_file = pysam.Samfile(params.perfect_reads_file, 'wb', header=bam_hdr)  # Write binary BAM with header

  bam_hdr = {'HD': {'VN': '1.4'},
             'SQ': [{'LN': ref_seq_len, 'SN': params.seq_header,
                     'SP': 'Simulated human, corrupted reads, by reads.py {:s}'.format(__version__)}]}
  corrupted_bam_file = pysam.Samfile(params.corrupted_reads_file, 'wb', header=bam_hdr)  # Write binary BAM with header


  prev_state = None
  read_count = 0
  while read_count < total_reads:
    corrupted_reads, perfect_reads, prev_state = \
      read_model.generate_reads(ref_seq, start=params.start, stop=params.stop,
                                num_reads=min(reads_per_block, total_reads - read_count),
                                prev_state=None, **params.args)

    save_reads_to_bam(corrupted_bam_file, corrupted_reads)
    save_reads_to_bam(perfect_bam_file, perfect_reads)
    read_count += len(corrupted_reads)   # corrupted_reads and perfect_reads should match exactly
    logger.debug('Write {:d} reads ({:d}%)'.format(read_count, int(100 * read_count / float(total_reads))))

  corrupted_bam_file.close()
  perfect_bam_file.close()
  f.close()
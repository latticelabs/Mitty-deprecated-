"""This module contains functions that generate simulated reads. The module can be called as a script as well. This is
useful for creating test data for MGR algorithms/data formats

Usage:
reads --ref=REF  [--start=START]  [--stop=STOP] [--coverage=COV] [--out=OUT] [-f] --paramfile=PFILE [--reads_per_block=BL] [-v]

Options:
  --ref=REF               The reference sequence in smalla format
  --start=START           Where to start taking reads (0,1) [default: 0.0]
  --stop=STOP             Where to stop taking reads (0,1) [default: 1.0]
  --coverage=COV          Coverage [default: 5]
  --out=OUT               Prefix of output file. [default: Data/simulated_reads]
                          The corrupted reads will be saved to simulated_reads.bam
                          The perfect reads will be saved to simulated_reads_perfect.bam
                          A text sidecar file simulated_reads.info will be saved with simulation parameters
  -f                      Write as FASTQ instead of BAM
  --paramfile=PFILE       Name for parameter file
  --reads_per_block=BL    Generate these many reads at a time (Adjust to machine resources). [default: 10000]
  -v                      Dump detailed logger messages

Notes:
1. The seq id of each read is the string 'rN:S1:S2' where N is the number of the read,
   S1 the start of the first read and S2 the start of the mate pair. Unpaired reads have no S2
2. The quality scores are in Phred scale (as specified in the SAM spec)

"""
#__version__ = '0.1.0'  # Fixed read lengths. Uniform coverage
#__version__ = '0.2.0'  # Paired end. Fixed read lengths. Uniform coverage
__version__ = '0.2.1'  # Plugin system implemented

import os
import imp
import json
import mmap
import docopt
import pysam  # Needed to write BAM files
import logging

logger = logging.getLogger(__name__)

# TODO: do we need seq_len?
def open_reads_files(out_prefix, seq_len, seq_header, bam=True):
  """Open output files (fastq or bam) and write headers as needed."""
  if bam:
    corrupted_out_fname = out_prefix + '.bam'
    perfect_out_fname = out_prefix + '_perfect.bam'
    bam_hdr = {'HD': {'VN': '1.4'},
               'SQ': [{'LN': seq_len, 'SN': seq_header,
                       'SP': 'Simulated perfect reads, by reads.py {:s}'.format(__version__)}]}
    perfect_fh = pysam.Samfile(perfect_out_fname, 'wb', header=bam_hdr)  # Write binary BAM with header
    bam_hdr['SQ'] = [{'LN': seq_len, 'SN': seq_header,
                     'SP': 'Simulated corrupted reads, by reads.py {:s}'.format(__version__)}]
    corrupted_fh = pysam.Samfile(corrupted_out_fname, 'wb', header=bam_hdr)  # Write binary BAM with header
  else:
    corrupted_out_fname = out_prefix + '.fastq'
    perfect_out_fname = out_prefix + '_perfect.fastq'
    #TODO: is there a fastq header?
    perfect_fh = open(perfect_out_fname, 'w')  # Write binary BAM with header
    corrupted_fh = open(corrupted_out_fname, 'w')  # Write binary BAM with header
  return perfect_fh, corrupted_fh


def save_reads(perfect_fh, corrupted_fh, perfect_reads, corrupted_reads, first_read_count, bam=True):
  """Wrapper around the functions that save the reads - calls appropriate functions depending on whether we want bam
   of fastq."""
  if bam:
    save_reads_to_bam(perfect_fh, perfect_reads, first_read_count)
    save_reads_to_bam(corrupted_fh, corrupted_reads, first_read_count)
  else:
    save_reads_to_fastq(perfect_fh, perfect_reads, first_read_count)
    save_reads_to_fastq(corrupted_fh, corrupted_reads, first_read_count)


def save_reads_to_bam(bam_file, reads, first_read_count):
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
    ar.qname = 'r{:d}:{:d}'.format(first_read_count+n, reads[n][0][2])
    ar.seq = reads[n][0][0]
    ar.qual = reads[n][0][1]
    if paired:
      ar.qname = ar.qname + ':{:d}'.format(reads[n][1][2])
      ar.flag = 0x41  # end1 0x01 flag has to be set to indicate multiple segments
    bam_file.write(ar)

    if paired:
      ar.seq = reads[n][1][0]
      ar.qual = reads[n][1][1]
      ar.flag = 0x81  # end2 0x01 flag has to be set to indicate multiple segments
      bam_file.write(ar)

def save_reads_to_fastq(fastq_file_handle, reads, first_read_count):
  """Given a list of sequences and their quality write the read data to an already opened text file. This saves data
  in interleaved form if paired reads are present

  Inputs:
    fastq_file_handle - file handles read for writing

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
    qname = 'r{:d}:{:d}'.format(first_read_count+n, reads[n][0][2])
    seq = reads[n][0][0]
    qual = reads[n][0][1]
    if paired:
      qname += ':{:d}'.format(reads[n][1][2])
    fastq_file_handle.write('@{:s}\n{:s}\n+\n{:s}\n'.format(qname, seq, qual))
    if paired:
      seq = reads[n][1][0]
      qual = reads[n][1][1]
      fastq_file_handle.write('@{:s}\n{:s}\n+\n{:s}\n'.format(qname, seq, qual))

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)
  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=logging.DEBUG)

  params = json.load(open(args['--paramfile'], 'r'))  #imp.load_source('params', args['--paramfile'], open(args['--paramfile'], 'r'))

  #Load the ref-seq smalla file
  f = open(args['--ref'], 'r+b')
  seq = mmap.mmap(f.fileno(), 0)
  seq_len = len(seq)
  seq_header = open(args['--ref'] + '.heada', 'r').readline()

  plugin_dir = os.path.join(os.path.dirname(__file__), 'Plugins', 'Reads')
  model_fname = os.path.join(plugin_dir, params['model'] + '_plugin.py')
  read_model = imp.load_source('readmodel', model_fname, open(model_fname, 'r'))

  start_reads = int(float(args['--start']) * seq_len)
  stop_reads = int(float(args['--stop']) * seq_len)
  subsequence_len = stop_reads - start_reads
  total_reads = int(float(args['--coverage']) * subsequence_len / float(read_model.average_read_len(**params['args'])))

  bam = not args['-f']  # bam is True if args['-f] is False
  perfect_fh, corrupted_fh = open_reads_files(args['--out'], seq_len, seq_header, bam)
  side_car_fname = args['--out'] + '.info'

  reads = read_model.read_generator(seq=seq,
                                    read_start=start_reads,
                                    read_stop=stop_reads,
                                    reads_per_call=int(args['--reads_per_block']),
                                    num_reads=total_reads,
                                    **params['args'])

  for corrupted_reads, perfect_reads, read_count in reads:
    save_reads(perfect_fh, corrupted_fh, perfect_reads, corrupted_reads, read_count - len(corrupted_reads), bam)
    logger.debug('Generated {:d} reads ({:d}%)'.format(read_count, int(100 * read_count / float(total_reads))))

  f.close()

  with open(side_car_fname,'w') as f:
    f.write('Command line\n-------------\n')
    f.write(json.dumps(args, indent=4))
    f.write('\n\nParameters\n------------\n')
    f.write(json.dumps(params, indent=4))

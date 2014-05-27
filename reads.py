"""This module contains functions that generate simulated reads. The module can be called as a script as well. This is
useful for creating test data for MGR algorithms/data formats

Usage:
reads [--paramfile=PFILE]  [--corrupt]  [--fastq] [--reads_per_block=BL]  [-v]

Options:
  --paramfile=PFILE       Name for parameter file (If none supplied, will read from stdin)
  --corrupt               Write out corrupted reads too.
  --fastq                 Write as FASTQ instead of BAM (simulated_reads.fastq)
  --reads_per_block=BL    Generate these many reads at a time (Adjust to machine resources). [default: 100000]
  -v                      Dump detailed logger messages

Notes:
1. The seq id of each read is the string 'rN:S1:S2' where N is the number of the read,
   S1 the start of the first read and S2 the start of the mate pair. Unpaired reads have no S2
2. The quality scores are in Phred scale (as specified in the SAM spec)
3. We supply the prefix of output file name in the parameter file . Say we set this as sim_reads.
   The perfect reads will be saved to sim_reads.bam (or sim_reads.fastq). If we ask for corrupted reads
   we will get the corrupted reads in the file sim_reads_c.fastq.
   A text sidecar file sim_reads.info will always be saved with simulation parameters.
"""

__explain__ = """
Example parameter file .json

{
    "input_sequences": ["mutated_1.smalla", "mutated_2.smalla"],
    "total_reads": [100, 100],
    "coverages": [2.0, 10.0],
    "is_this_ref_seq": false,
    "read_ranges": [[0.0, 1.0], [0.0, 1.0]],
    "output_file_prefix": "sim_reads",
    "read_model": "tiled_reads",
    "model_params": {
        "paired": false,
        "read_len": 100,
        "template_len": 250,
        "read_advance": 50
    }
}

Specify either total_reads OR coverage. If you specify both only the coverage will be honored
"""

__version__ = '0.3.0'

import os
import imp
import json
import mmap
import docopt
import numpy  # Needed for mmap of pos file
import pysam  # Needed to write BAM files
import logging

logger = logging.getLogger(__name__)


def open_reads_files(out_prefix, corrupted_reads=False, save_as_bam=True):
  """Open output files (fastq or bam) and write headers as needed."""
  file_handles = {}
  if save_as_bam:  # BAM
    perfect_reads_fname = out_prefix + '.bam'
    bam_hdr = {'HD': {'VN': '1.4'},
               'SQ': [{'LN': 1, 'SN': 'raw reads',
                       'SP': 'reads.py {:s}'.format(__version__)}]}
    file_handles['perfect'] = pysam.Samfile(perfect_reads_fname, 'wb', header=bam_hdr)  # Write binary BAM with header
    if corrupted_reads:
      corrupted_reads_fname = out_prefix + '_c.bam'
      bam_hdr['SQ'] = [{'LN': 1, 'SN': 'raw corrupted reads',
                        'SP': 'reads.py {:s}'.format(__version__)}]
      file_handles['corrupted'] = \
        pysam.Samfile(corrupted_reads_fname, 'wb', header=bam_hdr)  # Write binary BAM with header
  else:  # FASTQ
    perfect_reads_fname = out_prefix + '.fastq'
    file_handles['perfect'] = open(perfect_reads_fname, 'w')  # File handles for FASTQ
    if corrupted_reads:
      corrupted_reads_fname = out_prefix + '_c.fastq'
      file_handles['corrupted'] = open(corrupted_reads_fname, 'w')
  return file_handles


def roll_cigar(this_read, p_arr):
  """
  You can 'read along' with these tests from the Readme developers section

  Test a fully matching read
  >>> t_read = ['ATTG','~~~~', 0]; \
  p_arr = [1, 2, 3, 4, 5, 5, 5, 6, 9]; \
  roll_cigar(t_read, p_arr)
  (1, '4M')

  Test for read with insert
  >>> t_read = ['TTGT', '~~~~', 1]; \
  roll_cigar(t_read, p_arr)
  (2, '3M1I')

  Another test for read with insert
  >>> t_read = ['TGTT', '~~~~', 2]; \
  roll_cigar(t_read, p_arr)
  (3, '2M2I')

  Test for read with delete at end - should not show up in CIGAR
  >>> t_read = ['TTAC', '~~~~', 4]; \
  roll_cigar(t_read, p_arr)
  (5, '2I2M')

  Test for read spanning a deletion - should get a delete
  >>> t_read = ['ACAC', '~~~~', 0]; \
  p_arr = [1, 2, 5, 6, 7, 8, 9]; \
  roll_cigar(t_read, p_arr)
  (1, '2M2D2M')

  We actually missed this case: read with one matching base and then a delete
  >>> t_read = ['CACT', '~~~~', 1]; \
  roll_cigar(t_read, p_arr)
  (2, '1M2D3M')


  Test for an unmapped read: pos and cigars should be None
  >>> t_read = ['AATT', '~~~~', 2]; \
  p_arr = [1, 2, 3, 3, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9]; \
  roll_cigar(t_read, p_arr)
  (0, '')

  """
  mapped = False
  cigar = ''
  coord = this_read[2]
  counter = 0
  type = None
  for n in range(coord, coord + len(this_read[0])):
    dp = p_arr[n+1] - p_arr[n]
    if dp == 1:
      mapped = True  # As long as we have one M we are a mapped read
      if type != 'M':
        if counter > 0:  # Flush
          cigar += '{:d}{:s}'.format(counter, type)
          counter = 0
      type = 'M'
      counter += 1
    elif dp == 0:
      if type != 'I':
        if counter > 0:  # Flush
          cigar += '{:d}{:s}'.format(counter, type)
          counter = 0
      type = 'I'
      counter += 1
    elif dp > 1:
      mapped = True  # As long as we have one M we are a mapped read
      if type != 'M':
        if counter > 0:  # Flush
          cigar += '{:d}{:s}'.format(counter, type)
          counter = 0
      type = 'M'  # We need to set this because we could be at the start of a read and type = None still
      counter += 1
      cigar += '{:d}{:s}'.format(counter, type)
      type = 'D'
      counter = dp - 1

  if type != 'D':  # Flush all but 'D'. We only write D if we cross a D boundary
    cigar += '{:d}{:s}'.format(counter, type)

  if mapped:
    align_pos = p_arr[coord]
  else:
    align_pos = 0
    cigar = ''

  return align_pos, cigar


# The qname of an unpaired read is written as
# rN:POS:CIGAR
# while that of paired reads are written as
# rN:POS1:CIGAR1:POS2:CIGAR2
# This function fills out the POS and CIGAR in the qname
def roll(these_reads, pos_array):
  """Given a list of reads return us the reads with the coordinate replaced by a (POS,CIGAR) tuple. CIGAR, create the
  CIGAR, roll, get it? Oh, alright, I tried.
  Inputs:
                                 _________ ( seq_str, quality_str, coordinate)
    reads     -  [              /
                  [[ ... ], [ ... ]],
                  [[ ... ], [ ... ]], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads

    pos_array -  1 d array

  Output:
    No output - changes reads in place
  """
  def roll_null(these_reads):
    """Simple convenience function - for if we have null reads ."""
    paired = True if len(these_reads[0]) == 2 else False
    for this_read in these_reads:
      qname = '{:d}:{:d}M'.format(this_read[0][2] + 1, len(this_read[0][0]))
      if paired:
        qname += ':{:d}:{:d}M'.format(this_read[1][2] + 1, len(this_read[1][0]))
        this_read[1][2] = qname
      this_read[0][2] = qname


  if len(these_reads) == 0: return
  paired = True if len(these_reads[0]) == 2 else False
  # If these reads are from the reference sequence, then return CIGARs from the null read
  if pos_array is None:
    return roll_null(these_reads)

  for this_read in these_reads:
    pos, cigar = roll_cigar(this_read[0], pos_array)
    qname = '{:d}:{:s}'.format(pos, cigar)
    if paired:
      pos, cigar = roll_cigar(this_read[1], pos_array)
      qname += ':{:d}:{:s}'.format(pos, cigar)
      this_read[1][2] = qname
    this_read[0][2] = qname


def save_reads(fh, these_reads, offset, save_as_bam=True):
  if save_as_bam:
    save_reads_to_bam(fh, these_reads, offset)
  else:
    save_reads_to_fastq(fh, these_reads, offset)


# Didn't refactor these functions as that would involve another function call within a loop
#TODO: Implement option to write short qname
def save_reads_to_bam(bam_file, these_reads, first_read_serial):
  """
  Inputs:
                                 _________ ( seq_str, quality_str, qname)
    reads     -  [              /
                  [( ... ), ( ...)],
                  [( ... ), ( ...)], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads
  """
  if len(these_reads) == 0: return
  paired = True if len(these_reads[0]) == 2 else False
  for ser_no, this_read in enumerate(these_reads):
    ar = pysam.AlignedRead()
    ar.qname = 'r{:d}:'.format(first_read_serial + ser_no) + this_read[0][2]
    ar.seq = this_read[0][0]
    ar.qual = this_read[0][1]
    if paired:
      ar.flag = 0x41  # end1 0x01 flag has to be set to indicate multiple segments
    bam_file.write(ar)

    if paired:
      ar.seq = this_read[1][0]
      ar.qual = this_read[1][1]
      ar.flag = 0x81  # end2 0x01 flag has to be set to indicate multiple segments
      bam_file.write(ar)


def save_reads_to_fastq(fastq_file_handle, these_reads, first_read_serial):
  """Given a list of sequences and their quality write the read data to an already opened text file. This saves data
  in interleaved form if paired reads are present

  Inputs:
    fastq_file_handle - file handles read for writing

                                 _________ ( seq_str, quality_str, qname)
    reads     -  [              /
                  [( ... ), ( ...)],
                  [( ... ), ( ...)], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads

  """
  if len(these_reads) == 0: return
  paired = True if len(these_reads[0]) == 2 else False
  for ser_no, read in enumerate(these_reads):
    qname = 'r{:d}:'.format(first_read_serial + ser_no) + read[0][2]
    seq = read[0][0]
    qual = read[0][1]
    fastq_file_handle.write('@{:s}\n{:s}\n+\n{:s}\n'.format(qname, seq, qual))
    if paired:
      seq = read[1][0]
      qual = read[1][1]
      fastq_file_handle.write('@{:s}\n{:s}\n+\n{:s}\n'.format(qname, seq, qual))


def reads_from_a_sequence(seq_fname=None,
                          is_ref_seq=False,
                          read_count=100,
                          read_range=[0.0, 1.0],
                          read_model=None,
                          model_params=None,
                          read_file_handles=None,
                          write_corrupted=False,
                          save_as_bam=True,
                          reads_per_call=10000):

  with open(seq_fname, 'rb') as f_seq:
    seq = mmap.mmap(f_seq.fileno(), 0, access=mmap.ACCESS_READ)  # sequence the reads come from
    seq_len = len(seq)
    if is_ref_seq:
      seq_pos = None  # Only roll uses this and will assume reads from ref_seq if this is None
      logger.debug('Reference sequence, not searching for a pos file')
    else:
      seq_pos_filename = seq_fname + '.pos'
      if os.path.exists(seq_pos_filename):
        seq_pos = numpy.memmap(seq_pos_filename, dtype='int32', mode='r', shape=(seq_len + 1,))  # Our relevant pos file
      else:
        logger.error('Could not find POS file! Quitting')
        return

    start_reads = int(read_range[0] * seq_len)
    stop_reads = int(read_range[1] * seq_len)
    reads = read_model.read_generator(seq=seq,
                                      read_start=start_reads,
                                      read_stop=stop_reads,
                                      num_reads=read_count,
                                      reads_per_call=reads_per_call,
                                      **model_params)
    current_read_count = 0
    for perfect_reads in reads:
      roll(perfect_reads, seq_pos)  # The function modifies perfect_reads in place to fill out the POS and CIGAR
      save_reads(read_file_handles['perfect'], perfect_reads, current_read_count, save_as_bam)
      #save_reads needs current_read_count because we put in a read serial number in the qname
      if write_corrupted:
        save_reads(read_file_handles['corrupted'], read_model.corrupt_reads(perfect_reads, **model_params),
                   current_read_count, save_as_bam)
      current_read_count += len(perfect_reads)
      logger.debug('Generated {:d} reads ({:d}%)'.
                   format(current_read_count, int(100 * current_read_count / float(read_count))))


def main(args):
  # Read parameter file from stdin if no name is given
  param_string = sys.stdin.read() if args['--paramfile'] is None else open(args['--paramfile'], 'r').read()
  params = json.loads(param_string)

  # Load the read model from the plugins directory
  plugin_dir = os.path.join(os.path.dirname(__file__), 'Plugins', 'Reads')
  model_fname = os.path.join(plugin_dir, params['read_model'] + '_plugin.py')
  read_model = imp.load_source('readmodel', model_fname, open(model_fname, 'r'))

  # Generate a dictionary of file handles for perfect and corrupted reads (if needed)
  save_as_bam = not args['--fastq']  # bam is True if args['-f] is False
  write_corrupted = args['--corrupt']  # If True, corrupted reads will be written out
  read_file_handles = open_reads_files(params['output_file_prefix'], write_corrupted, save_as_bam)

  #if coverage is specified fill out total_reads based on sequence lengths etc.
  if 'coverages' in params:
    params['total_reads'] = \
      [int(os.stat(seq_name).st_size * float(cov) / read_model.average_read_len(params['model_params']))
       for seq_name, cov in zip(params['input_sequences'], params['coverages'])]

  # For each sequence in our input list generate reads and flush to file
  for seq_name, read_count, read_range in zip(params['input_sequences'], params['total_reads'], params['read_ranges']):
    reads_from_a_sequence(seq_fname=seq_name,
                          is_ref_seq=params['is_this_ref_seq'],
                          read_count=read_count,
                          read_range=read_range,
                          read_model=read_model,
                          model_params=params['model_params'],
                          read_file_handles=read_file_handles,
                          write_corrupted=write_corrupted,
                          save_as_bam=save_as_bam,
                          reads_per_call=int(args['--reads_per_block']))

  with open(params['output_file_prefix'] + '.info', 'w') as f:
    f.write('Command line\n-------------\n')
    f.write(json.dumps(args, indent=4))
    f.write('\n\nParameters\n------------\n')
    f.write(json.dumps(params, indent=4))


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  elif docopt.sys.argv[1] == 'test':
    import doctest
    doctest.testmod()
    exit(0)
  elif docopt.sys.argv[1] == 'explain':
    print __explain__
    exit(0)
  else:
    arguments = docopt.docopt(__doc__, version=__version__)
  level = logging.DEBUG if arguments['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  main(arguments)

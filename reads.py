"""This module contains functions that generate simulated reads. The module can be called as a script as well. This is
useful for creating test data for MGR algorithms/data formats

Usage:
reads --ref=REF --paramfile=PFILE [--start=START]  [--stop=STOP] [--coverage=COV] [--out=OUT] [-c] [-f] [--shortqname] [--reads_per_block=BL] [-v]

Options:
  --ref=REF               The reference sequence in smalla format
  --start=START           Where to start taking reads (0,1) [default: 0.0]
  --stop=STOP             Where to stop taking reads (0,1) [default: 1.0]
  --coverage=COV          Coverage [default: 5]
  --out=OUT               Prefix of output file. [default: Data/simulated_reads]
                          The perfect reads will be saved to simulated_reads.bam
                          A text sidecar file simulated_reads.info will be saved with simulation parameters
  -c                      Write corrupted reads (simulated_reads_c.bam)
  -f                      Write as FASTQ instead of BAM (simulated_reads.fastq)
  --shortqname            Instead of writing the POS and CIGAR into the qname, reads.py will only write the read serial
                          number. Some tools, such as Tablet and IGV can be crashed if len(qname) > 254
  --paramfile=PFILE       Name for parameter file
  --reads_per_block=BL    Generate these many reads at a time (Adjust to machine resources). [default: 10000]
  -v                      Dump detailed logger messages

Notes:
1. The seq id of each read is the string 'rN:S1:S2' where N is the number of the read,
   S1 the start of the first read and S2 the start of the mate pair. Unpaired reads have no S2
2. The quality scores are in Phred scale (as specified in the SAM spec)

"""
__version__ = '0.3.0'

import os
import imp
import json
import mmap
import docopt
import pysam  # Needed to write BAM files
import logging

logger = logging.getLogger(__name__)


def open_reads_files(out_prefix, seq_len, seq_header, corrupted_reads=False, save_as_bam=True):
  """Open output files (fastq or bam) and write headers as needed."""
  file_handles = {}
  if save_as_bam:  # BAM
    perfect_reads_fname = out_prefix + '.bam'
    bam_hdr = {'HD': {'VN': '1.4'},
               'SQ': [{'LN': seq_len, 'SN': seq_header,
                       'SP': 'Simulated perfect reads, by reads.py {:s}'.format(__version__)}]}
    file_handles['perfect'] = pysam.Samfile(perfect_reads_fname, 'wb', header=bam_hdr)  # Write binary BAM with header
    if corrupted_reads:
      corrupted_reads_fname = out_prefix + '_c.bam'
      bam_hdr['SQ'] = [{'LN': seq_len, 'SN': seq_header,
                        'SP': 'Simulated corrupted reads, by reads.py {:s}'.format(__version__)}]
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
  You can 'read along' with this tests from the Readme

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
      mapped = True  # As long as we have one I we are a mapped read
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
      mapped = True  # As long as we have one I we are a mapped read
      if type != 'M':
        if counter > 0:  # Flush
          cigar += '{:d}{:s}'.format(counter, type)
          counter = 0
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

    reads    -  [
                  [(pos, cigar), (pos, cigar)], -> inner list = 2 elements if paired reads, 1 otherwise
                       .
                       .
                       .
                 ] -> outer list = number of reads

  """
  def roll_null(these_reads):
    """Simple convenience function - for if we have null reads ."""
    paired = True if len(these_reads[0]) == 2 else False
    for this_read in these_reads:
      qname = '{:d}:{:d}M'.format(this_read[0][2], len(this_read[0][0]))
      if paired:
        qname += ':{:d}:{:d}M'.format(this_read[1][2], len(this_read[1][0]))
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


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  elif docopt.sys.argv[1] == 'test':
    import sys
    import doctest
    doctest.testmod()
    sys.exit()
  else:
    args = docopt.docopt(__doc__, version=__version__)
  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  params = json.load(open(args['--paramfile'], 'r'))  # The parameter file in json format

  # Load the read model from the plugins directory
  plugin_dir = os.path.join(os.path.dirname(__file__), 'Plugins', 'Reads')
  model_fname = os.path.join(plugin_dir, params['model'] + '_plugin.py')
  read_model = imp.load_source('readmodel', model_fname, open(model_fname, 'r'))

  #Load the ref-seq smalla file and read the header side car
  with open(args['--ref'], 'r+b') as f_ref:
    seq = mmap.mmap(f_ref.fileno(), 0)  # Our reference sequence
    seq_len = len(seq)
    seq_header = open(args['--ref'] + '.heada', 'r').readline()

    ref_pos_fname = args['--ref'] + '.pos'
    if os.path.exists(ref_pos_fname):
      f_ref_pos = open(ref_pos_fname, 'r+b')
      seq_pos = mmap.mmap(f_ref_pos.fileno(), 0)  # Our relevant pos file
    else:
      seq_pos = None  # Only roll uses this and will assume reads from ref_seq if this is None

    start_reads = int(float(args['--start']) * seq_len)
    stop_reads = int(float(args['--stop']) * seq_len)
    total_reads_required = int(float(args['--coverage']) * (stop_reads - start_reads) /
                               float(read_model.average_read_len(**params['args'])))

    # Generate a dictionary of file handles for perfect and corrupted reads (if needed)
    bam = not args['-f']  # bam is True if args['-f] is False
    write_corrupted = args['-c']  # If True, corrupted reads will be written out
    read_file_handles = open_reads_files(args['--out'], seq_len, seq_header, write_corrupted, bam)

    reads = read_model.read_generator(seq=seq,
                                      read_start=start_reads,
                                      read_stop=stop_reads,
                                      reads_per_call=int(args['--reads_per_block']),
                                      num_reads=total_reads_required,
                                      **params['args'])
    current_read_count = 0
    for perfect_reads in reads:
      roll(perfect_reads, seq_pos)  # The function modifies perfect_reads in place to fill out the POS and CIGAR
      save_reads(read_file_handles['perfect'], perfect_reads, current_read_count, bam)
      #save_reads needs current_read_count because we put in a read serial number in the qname
      if write_corrupted:
        save_reads(read_file_handles['corrupted'], read_model.corrupt_reads(perfect_reads, **params['args']),
                   current_read_count, bam)
      current_read_count += len(perfect_reads)
      logger.debug('Generated {:d} reads ({:d}%)'.
                   format(current_read_count, int(100 * current_read_count / float(total_reads_required))))

  with open(args['--out'] + '.info','w') as f:
    f.write('Command line\n-------------\n')
    f.write(json.dumps(args, indent=4))
    f.write('\n\nParameters\n------------\n')
    f.write(json.dumps(params, indent=4))

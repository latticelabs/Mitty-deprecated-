"""Functions for simplified but fast IO of FASTA formatted files.

TODO: Unittest or some other testing framework. Main problem: how to handle data - make small data sets and save with
the code - download small genome fasta file from NCBI

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""


def int_read_fasta(file_handle):
  """This is a simplified FASTA reader that runs slowa than biopython but converts sequences into ints.
  This expects the fasta file to carry only one sequence.
  """
  import numpy
  file_handle.seek(0,2)
  file_size = file_handle.tell()
  file_handle.seek(0)
  seq = numpy.empty(file_size, dtype='int8')
  header = file_handle.readline()[:-1]
  index = 0
  for line in file_handle.readlines():
    seq_len = len(line) - 1
    seq[index:index+seq_len] = numpy.frombuffer(line[:-1],dtype='int8')
    index += seq_len
  return header, seq[:index]


def fast_read_fasta(file_handle):
  """This is a simplified FASTA reader that runs fasta than biopython. This expects the fasta file to carry only one
  sequence.
  Inputs:
    file_handle       - an open file handle

  >>> header, seq = fast_read_fasta(open('porcine_circovirus.fa','r'))
  >>> header
  'gi|52547303|gb|AY735451.1| Porcine circovirus isolate Hebei capsid protein gene, complete cds'
  >>> seq[-20:]
  bytearray(b'ACCCCCCACTTAACCCTTAG')
  """
  seq = bytearray()
  header = file_handle.readline()[1:-1]  # Get rid of > sign
  for line in file_handle.readlines():
    seq += line.rstrip()
  return header, seq

def fast_write_fasta(file_handle, seq, header=None, width=71):
  """Write out single sequence with header with a appropriate line breaks.
  Inputs:
    file_handle       - an open file handle
    seq               - the sequence
    header            - if None then no header will be written. Good for batch writing long sequences
    width             - when printing the sequence limit this width
  """
  if header is not None: file_handle.write('>' + header + '\n')
  for index in range(0, len(seq), width):
    file_handle.write(seq[index:index+width] + '\n')


def time_test():
  """Runs the fast_read_fasta function to compare it against BioPython."""
  import timeit

  tm = timeit.Timer(stmt="seq = SeqIO.read('porcine_circovirus.fa', 'fasta')", setup='from Bio import SeqIO')
  print 'BioPython'
  print tm.timeit(number=2000)

  tm = timeit.Timer(stmt="header, seq = seqio.fast_read_fasta(open('porcine_circovirus.fa'))", setup='import seqio')
  print 'fast_read_fasta'
  print tm.timeit(number=2000)


def write_test():
  """Reads and then writes out a fasta file and then uses BioPython to spot compare the original and written sequences
  to test if they are same.

  >>> write_test()
  ACCCCCCACTTAACCCTTAG

  """
  from Bio import SeqIO
  header, seq = fast_read_fasta(open('porcine_circovirus.fa'))
  fast_write_fasta(open('test.fa','w'), seq, header, width=70)
  seq_orig = SeqIO.read('porcine_circovirus.fa','fasta')
  seq_written = SeqIO.read('test.fa','fasta')
  #print seq_orig.seq[-20:]
  print seq_written.seq[-20:]


if __name__ == "__main__":
  time_test()
  import doctest
  doctest.testmod()
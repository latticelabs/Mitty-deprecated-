"""Script to convert from fasta to smalla files

Usage:
converta <fasta> <smalla> [--block_size=BS]

Options:
  fasta              Input file name
  smalla             Output file name
  --block_size=BS    Block size [default: 1000000]

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
import docopt
__version__ = '0.1.0'

def get_fasta_header(fname):
  """Read FASTA header

  >>> print get_fasta_header(fname='Data/porcine_circovirus.fa')
  gi|52547303|gb|AY735451.1| Porcine circovirus isolate Hebei capsid protein gene, complete cds
  """
  line = open(fname, 'r').readline().strip()
  if len(line): line = line[1:]  # get rid of the leading >
  return line


def block_read_fasta(fname, block_size=1000):
  """Read FASTA file in blocks. For now this simply discards any lines begining with '>'

  >>> print ''.join([s for s in block_read_fasta(fname='Data/porcine_circovirus.fa', block_size=1)])
  ATGACGTATCCAAGGAGGCGTTACCGGAGAAGAAGACACCGCCCCCGCAGCCATCTTGGCCAGATCCTCCGCCGCCGCCCCTGGCTCGTCCACCCCCGCCACCGTTACCGCTGGAGAAGGAAAAACGGCATCTTCAACACCCGCCTCTCCCGCACCTTCGGATATACTATCAAGCGAACCACAGTCAAAACGCCCTCCTGGGCGGTGGACATGATGAGATTCAATATTAATGACTTTCTTCCCCCAGGAGGGGGCTCAAACCCCCGCTCTGTGCCCTTTGAATACTACAGAATAAGAAAGGTTAAGGTTGAATTCTGGCCCTGCTCCCCGATCACCCAGGGTGACAGGGGAGTGGGCTCCAGTGCTGTTATTCTAGATGATAACTTTGTAACAAAGGCCACAGCCCTCACCTATGACCCCTATGTAAACTACTCCTCCCGCCATACCATAACCCAGCCCTTCTCCTACCACTCCCGCTACTTTACCCCCAAACCTGTCCTAGATTCCACTATTGATTACTTCCAACCAAACAACAAAAGAAATCAGCTGTGGCTGAGACTACAAACTGCTGGAAATGTAGACCACGTAGGCCTCGGCACTGCGTTCGAAAACAGTATATACGACCAGGAATACAATATCCGTGTAACCATGTATGTACAATTCAGAGAATTTAATCTTAAAGACCCCCCACTTAACCCTTAG
  >>> print ''.join([s for s in block_read_fasta(fname='Data/porcine_circovirus.fa', block_size=13)])
  ATGACGTATCCAAGGAGGCGTTACCGGAGAAGAAGACACCGCCCCCGCAGCCATCTTGGCCAGATCCTCCGCCGCCGCCCCTGGCTCGTCCACCCCCGCCACCGTTACCGCTGGAGAAGGAAAAACGGCATCTTCAACACCCGCCTCTCCCGCACCTTCGGATATACTATCAAGCGAACCACAGTCAAAACGCCCTCCTGGGCGGTGGACATGATGAGATTCAATATTAATGACTTTCTTCCCCCAGGAGGGGGCTCAAACCCCCGCTCTGTGCCCTTTGAATACTACAGAATAAGAAAGGTTAAGGTTGAATTCTGGCCCTGCTCCCCGATCACCCAGGGTGACAGGGGAGTGGGCTCCAGTGCTGTTATTCTAGATGATAACTTTGTAACAAAGGCCACAGCCCTCACCTATGACCCCTATGTAAACTACTCCTCCCGCCATACCATAACCCAGCCCTTCTCCTACCACTCCCGCTACTTTACCCCCAAACCTGTCCTAGATTCCACTATTGATTACTTCCAACCAAACAACAAAAGAAATCAGCTGTGGCTGAGACTACAAACTGCTGGAAATGTAGACCACGTAGGCCTCGGCACTGCGTTCGAAAACAGTATATACGACCAGGAATACAATATCCGTGTAACCATGTATGTACAATTCAGAGAATTTAATCTTAAAGACCCCCCACTTAACCCTTAG
  >>> print ''.join([s for s in block_read_fasta(fname='Data/porcine_circovirus.fa', block_size=1000)])
  ATGACGTATCCAAGGAGGCGTTACCGGAGAAGAAGACACCGCCCCCGCAGCCATCTTGGCCAGATCCTCCGCCGCCGCCCCTGGCTCGTCCACCCCCGCCACCGTTACCGCTGGAGAAGGAAAAACGGCATCTTCAACACCCGCCTCTCCCGCACCTTCGGATATACTATCAAGCGAACCACAGTCAAAACGCCCTCCTGGGCGGTGGACATGATGAGATTCAATATTAATGACTTTCTTCCCCCAGGAGGGGGCTCAAACCCCCGCTCTGTGCCCTTTGAATACTACAGAATAAGAAAGGTTAAGGTTGAATTCTGGCCCTGCTCCCCGATCACCCAGGGTGACAGGGGAGTGGGCTCCAGTGCTGTTATTCTAGATGATAACTTTGTAACAAAGGCCACAGCCCTCACCTATGACCCCTATGTAAACTACTCCTCCCGCCATACCATAACCCAGCCCTTCTCCTACCACTCCCGCTACTTTACCCCCAAACCTGTCCTAGATTCCACTATTGATTACTTCCAACCAAACAACAAAAGAAATCAGCTGTGGCTGAGACTACAAACTGCTGGAAATGTAGACCACGTAGGCCTCGGCACTGCGTTCGAAAACAGTATATACGACCAGGAATACAATATCCGTGTAACCATGTATGTACAATTCAGAGAATTTAATCTTAAAGACCCCCCACTTAACCCTTAG
  """
  seq_buff = ''
  seq_buff_len = 0
  with open(fname, 'r') as f:
    for line in f.readlines():
      if line[0] == '>':

        continue  # A seq id line
      line = line.strip()
      seq_buff_len += len(line)
      seq_buff += line
      while seq_buff_len > block_size:  # Time to spit it out
        seq = seq_buff[:block_size]
        seq_buff = seq_buff[block_size:]
        seq_buff_len -= block_size
        yield seq

    if seq_buff_len:
      yield seq_buff  # The last tail after the file ends


def fasta_to_smalla(fasta_fname, smalla_fname, block_size=1000):
  with open(smalla_fname + '.heada', 'w') as f:
    f.write(get_fasta_header(args['<fasta>']))

  with open(smalla_fname, 'w') as f:
    for s in block_read_fasta(fname=fasta_fname, block_size=block_size):
      f.write(s)


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
  fasta_to_smalla(args['<fasta>'], args['<smalla>'], int(args['--block_size']))

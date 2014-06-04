"""
For efficiency purposes we strip the header and all new lines out of the original .fasta file. The resulting file is
called a .smalla file and is what the rest of the tools use. The header is saved into a .heada file and carries the
sequence header on one line and the sequence length on the other.

Usage:
converta  <fasta>  <smalla_prefix>
converta test [-v]

Options:
  fasta              Input file name
  smalla_prefix      Output file prefix

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
import docopt
__version__ = '0.1.0'


def fasta_to_smalla(fasta_fname, smalla_prefix):
  """Flush each sequence in the file into a separate smalla file, write headers into side-car files
  >>> import tempfile, os
  >>> seq1 = 'ATGACGTATCCAAGGAGGCGTTACCGGAGAAGAAGACACCGCCCCC'
  >>> seq2 = 'GCAGCCATCTTGGCCAGATCCTCCGCCGCCGCCCCTGGCTCGTC'
  >>> seq3 = 'ATGACGTATCCAAGGAGGCCCTCCGCCGCCGCCCCTGGCTCGTCCA'
  >>> tdir = tempfile.gettempdir()
  >>> open(os.path.join(tdir, 'test.fa'), 'w').write('>hdr1\\nATGACGTATCCAA\\nGGAGGCGTTACCGGAGAA\\nGAAG\\nA\\nCACCGCCCCC\\n>hdr2\\nGCAGCCA\\nTC\\nT\\nTGGCCAGATCCTCCGCCGCCGCCCC\\nTGGCTCGTC>hdr3\\nATGACGTATCCAAGGAGGCCCTCCGCCGCCGCCCCTGGCTCGTCCA')
  >>> fasta_to_smalla(os.path.join(tdir, 'test.fa'), os.path.join(tdir,'test'))
  >>> open(os.path.join(tdir, 'test_0.smalla.heada'), 'r').read()
  'hdr1\\n46'
  >>> open(os.path.join(tdir, 'test_1.smalla.heada'), 'r').read()
  'hdr2\\n44'
  >>> open(os.path.join(tdir, 'test_2.smalla.heada'), 'r').read()
  'hdr3\\n46'
  >>> open(os.path.join(tdir, 'test_0.smalla'), 'r').read() == seq1
  True
  >>> open(os.path.join(tdir, 'test_1.smalla'), 'r').read() == seq2
  True
  >>> open(os.path.join(tdir, 'test_2.smalla'), 'r').read() == seq3
  True
  """
  with open(fasta_fname, 'r') as fin:
    read_a_byte = fin.read
    seq_count = 0  # How many sequences are we at now?
    fout = None
    while True:
      byte = read_a_byte(1)
      if not byte or byte == '>':  # We've gotten to a new sequence header or the end of the file
                                   # time to flush the prev sequence and start a new one
        if fout:  # We were reading a sequence
          with open(smalla_fname + '.heada', 'w') as f:  # Write the header (.smalla.heada file)
            f.write('{:s}\n{:d}'.format(seq_name, fout.tell()))  # Header name, size of the written sequence
          fout.close()  # Close the sequence (.smalla) file
        if not byte:
          break
        smalla_fname = '{:s}_{:d}.smalla'.format(smalla_prefix, seq_count)
        seq_name = fin.readline()[:-1]
        fout = open(smalla_fname, 'w')
        seq_count += 1
      elif byte != '\n':  # Part of the sequence
        fout.write(byte)

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)
  if args['test']:
    import doctest
    doctest.testmod()
    exit(0)

  fasta_to_smalla(args['<fasta>'], args['<smalla_prefix>'])

"""
For efficiency purposes we strip the header and all new lines out of the original .fasta file. The resulting file is
called a .smalla file and is what the rest of the tools use. The header is saved into a .heada file and carries the
sequence header on one line and the sequence length on the other.

Usage:
converta  --datadir=DDIR   --descfile=DSCF   --smalla=SMLA  [-v]
converta tofasta --smalla=SMLA  --fasta=FASTA [-v]
converta explain
converta test [-v]

Options:
  --datadir=DDIR     The data root for all the input fasta files
  --descfile=DSCF    Description file in .json format
  --smalla=SMLA      Name of smalla file to save to
  tofasta            Take a smalla file and create a combined fasta file from all the sequences
  --fasta=FASTA      Output fasta file name
  -v                 Be verbose when you do things
  explain            Print out formats of the description file and the smalla file
  test               Run tests

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
import json
import struct
import os
import docopt
__version__ = '0.2.0'

__explain__ = """
Example description file

{
    "header": {
        "species": "Test chimera",
        "chromosome count": 2,
        "ploidy of data": 2
    },
    "chromosome names": {
        "1": "chr1",
        "2": "chr2"
    },
    "files": {
        "porcine_circovirus.fa": [
            {
                "chromosome number": 1,
                "chromosome copy": 1
            }
        ],
        "adenovirus.fa": [
            {
                "chromosome number": 2,
                "chromosome copy": 2
            }
        ],
        "chimera.fa": [
            {
                "chromosome number": 1,
                "chromosome copy": 2
            },
            {
                "chromosome number": 2,
                "chromosome copy": 1
            }
        ]
    }
}

Note that though converta can handle such out of order descriptions, it is better for readability and bug checking your
data to keep the numbering in order

.smalla format

    Header-----------------------------------------------
    [char10]   - version string of smalla format
    [char255]  - Human readable species name (string)
    [uint16]   - number of chromosomes max 65535
    [uint8]    - ploidy of the *data* (1,2,3 ...) max 255

    Index-------------------------------------------------
    For each chromosome, each copy --------------------------
        [uint16]  - chromosome number (1,2,3,4 ...)
        [uint16]  - chromosome copy (1,2,3 ...)
        [char255] - Human readable description
        [char255] - NCBI accession string (if any)
        [uint32]  - start byte of data in this file
        [uint32]  - length of sequence

    Data--------------------------------------------------
    For each chromosome copy --------------------------
      [uchar]   - Nucleotide data
        ...
"""

header_fmt = '10s 255s H B'
index_fmt = 'H H 255s 255s L L'


def create_smalla(desc_file, data_dir, smalla_file):
  with open(desc_file, 'r') as f:
    desc = json.load(f)

  with open(smalla_file, 'w') as fp:
    write_header(fp, desc)
    seq_count = 0

    for fasta_file, seq_info in desc['files'].iteritems():
      pass






def write_header(fp, desc):
  """Given a freshly opened file, write the file header from the information given in the description file."""
  fp.seek(0)
  values = (__version__, desc['header']['species'], desc['header']['chromosome count'], desc['header']['ploidy of data'])
  fp.write(struct.pack(header_fmt, *values))


def write_index(fp, chrom_list):
  """We do this only at the end, when we have all the accession strings, location and lengths of sequences and so have
  built up a complete index"""
  fp.seek(struct.calcsize(header_fmt))
  for chrom in chrom_list:
    fp.write(struct.pack(index_fmt, *(chrom['chromosome number'], chrom['copy'], chrom['chromosome description'],
                                      chrom['id string'],
                                      chrom['start byte of sequence data'], chrom['length of sequence'])))




#def add_sequence(seq_header, )


def build_smalla(fp, desc, data_dir):
  """Given file with ."""




def copy_fasta(fasta_fp, smalla_fp):
  """Given file pointers to a fasta file and a smalla file copy one sequence over.
  Inputs:
    fasta_fp       - fasta file pointer
    smalla_fp      - our smalla file, pointing at where we should be writing sequence information to

  >>> import io

  Testing basic reading of one sequence
  >>> with io.BytesIO() as smalla_fp, open('Data/chimera.fa','r') as fasta_fp:
  ...     print copy_fasta(fasta_fp, smalla_fp)
  ...     _ = smalla_fp.seek(0)
  ...     print smalla_fp.read(70)
  ('gi|448872318|gb|KC140233.1| Human herpesvirus 1 isolate 0215T/GS/CHN/2008 glycoprotein G (US4) gene, complete cds', 717L)
  ATGTCGCCGGGCGCCATGCGTGCCGTTGTTCCCATTATCCCATTCCTTTTGGTTCTTGTCGGTGTATCGG

  Testing reading of all sequences in a file
  >>> with io.BytesIO() as smalla_fp, open('Data/chimera.fa','r') as fasta_fp:
  ...     seq_l = 1
  ...     while seq_l:
  ...         seq_h, seq_l = copy_fasta(fasta_fp, smalla_fp)
  ...         print seq_h, seq_l
  ...     _ = smalla_fp.seek(717)
  ...     print smalla_fp.read(70)
  gi|448872318|gb|KC140233.1| Human herpesvirus 1 isolate 0215T/GS/CHN/2008 glycoprotein G (US4) gene, complete cds 717
  gi|9626078|ref|NC_001358.1| Parvovirus H1, complete genome 5176
  None 0
  CATTTTTAGAACTGACCAACCATGTTCACGCAAGTGACGTGATGACGCGCGCTGCGCGCGCTGCCTTCGG
  """
  read_a_byte = fasta_fp.read
  seq_name = None
  seq_start_pos = None
  byte = read_a_byte(1)
  while byte:
    if byte == '>':  # We've gotten to a new sequence header
      if seq_name:  # We have been copying over a sequence, so we need to report back to caller
        fasta_fp.seek(-1, 1)  # Pretend this never happened
        break
      else:  # We are seeing a sequence for the first time this call
        seq_name = fasta_fp.readline()[:-1]
        seq_start_pos = smalla_fp.tell()
    elif byte != '\n':  # Part of the sequence, copy it to the smalla file
      smalla_fp.write(byte)
    byte = read_a_byte(1)
  seq_len = smalla_fp.tell() - seq_start_pos if seq_start_pos is not None else 0
  return seq_name, seq_len




if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)
  if args['explain']:
    print __explain__
    exit(0)
  if args['test']:
    import doctest
    doctest.testmod()
    exit(0)

  fasta_to_smalla(args['<fasta>'], args['<smalla_prefix>'])

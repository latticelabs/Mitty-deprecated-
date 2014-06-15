"""
This module defines a whole genome file (.wg) that allows us to store multiple sequences in an organized fashion. It
also defines and interface to
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
import numpy
import docopt
__version__ = '0.2.0'

__explain__ = """
Example description file

{
    "header": {
        "species": "Test Chimera",
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


def read_header(smalla_fp):
  """This function preserves the file pointer location."""
  cp = smalla_fp.tell()
  smalla_fp.seek(0)
  values = struct.unpack(header_fmt, smalla_fp.read(struct.calcsize(header_fmt)))
  smalla_fp.seek(cp)
  return {
    'smalla version': values[0],
    'species': values[1],
    'chromosome count': values[2],
    'ploidy': values[3]
  }


def read_index(smalla_fp, header):
  """This function preserves the file pointer location."""
  def read_index_entry(fp):
    values = struct.unpack(index_fmt, fp.read(struct.calcsize(index_fmt)))
    return {
      'chromosome number': values[0],
      'chromosome copy': values[1],
      'chromosome description': values[2].strip('\0'),
      'sequence id': values[3].strip('\0'),
      'start byte of sequence data': values[4],
      'length of sequence': values[5]
    }
  cp = smalla_fp.tell()
  smalla_fp.seek(struct.calcsize(header_fmt))
  index = {}
  for n in range(header['chromosome count'] * header['ploidy']):
    this_index = read_index_entry(smalla_fp)
    index['{:d}:{:d}'.format(this_index['chromosome number'], this_index['chromosome copy'])] = this_index
  smalla_fp.seek(cp)
  return index


def get_sequence(smalla_fp, index, chrom, copy):
  offset = index['{:d}:{:d}'.format(chrom, copy)]['start byte of sequence data']
  seq_len = index['{:d}:{:d}'.format(chrom, copy)]['length of sequence']
  return numpy.memmap(smalla_fp, dtype='c', mode='r', offset=offset, shape=(seq_len,))


def create_smalla(desc, data_dir, smalla_fp):
  """Main function to create a smalla file from a desc file and directory of fasta files.

  >>> import io, tempfile
  >>> with open('Data/wg_desc.json', 'r') as f:
  ...    desc = json.load(f)
  >>> smalla_fp = tempfile.TemporaryFile()
  >>> create_smalla(desc, 'Data/', smalla_fp)
  >>> hdr = read_header(smalla_fp)
  >>> print hdr['species']
  Test Chimera
  >>> print hdr['ploidy']
  2
  >>> index = read_index(smalla_fp, hdr)
  >>> print index['1:2']['sequence id']
  gi|448872318|gb|KC140233.1| Human herpesvirus 1 isolate 0215T/GS/CHN/2008 glycoprotein G (US4) gene, complete cds
  >>> print index['1:1']['sequence id']
  gi|52547303|gb|AY735451.1| Porcine circovirus isolate Hebei capsid protein gene, complete cds
  >>> seq = get_sequence(smalla_fp, index, 1, 1)
  >>> seq[:10]
  ATGACGTATCCAAGGAGGCGTTACCGGAGAAGAAGACACCGCCCCCGCAGCCATCTTGGCCAGATCCTCCGCCGCCGCCCCTGGCTCGTCCACCCCCGCCACCGTTACCGCTGGAGAAGGAAAAACGGCATCTTCAACACCCGCCTCTCCCGCACCTTCGGATATACTATCAAGCGAACCACAGTCAAAACGCCCTCCTGGGCGGTGGACATGATGAGATTCAATATTAATGACTTTCTTCCCCCAGGAGGGGGCTCAAACCCCCGCTCTGTGCCCTTTGAATACTACAGAATAAGAAAGGTTAAGGTTGAATTCTGGCCCTGCTCCCCGATCACCCAGGGTGACAGGGGAGTGGGCTCCAGTGCTGTTATTCTAGATGATAACTTTGTAACAAAGGCCACAGCCCTCACCTATGACCCCTATGTAAACTACTCCTCCCGCCATACCATAACCCAGCCCTTCTCCTACCACTCCCGCTACTTTACCCCCAAACCTGTCCTAGATTCCACTATTGATTACTTCCAACCAAACAACAAAAGAAATCAGCTGTGGCTGAGACTACAAACTGCTGGAAATGTAGACCACGTAGGCCTCGGCACTGCGTTCGAAAACAGTATATACGACCAGGAATACAATATCCGTGTAACCATGTATGTACAATTCAGAGAATTTAATCTTAAAGACCCCCCACTTAACCCTTAG
  """
  write_header(smalla_fp, desc)
  pointer_to_index = struct.calcsize(header_fmt)  # The index starts right after the header
  # The sequence data starts after the index. There are chromosome count x ploidy of data entries in the index
  pointer_to_sequence = \
    struct.calcsize(index_fmt) * desc['header']['chromosome count'] * desc['header']['ploidy of data'] \
    + struct.calcsize(header_fmt)
  for fasta_file, seq_info in desc['files'].iteritems():
    with open(os.path.join(data_dir, fasta_file), 'r') as fasta_fp:
      seq_cnt = 0
      smalla_fp.seek(pointer_to_sequence)
      seq_h, seq_l = copy_fasta(fasta_fp, smalla_fp)
      while seq_l:
        seq_info[seq_cnt]['start byte of sequence data'] = pointer_to_sequence  # We originally started to write from here
        seq_info[seq_cnt]['chromosome description'] = desc['chromosome names'][str(seq_info[seq_cnt]['chromosome number'])]
        seq_info[seq_cnt]['sequence id'] = seq_h
        seq_info[seq_cnt]['length of sequence'] = seq_l
        pointer_to_sequence = smalla_fp.tell()  # Next sequence will start from here
        smalla_fp.seek(pointer_to_index)
        append_index(smalla_fp, seq_info[seq_cnt])
        pointer_to_index = smalla_fp.tell()
        smalla_fp.seek(pointer_to_sequence)
        seq_h, seq_l = copy_fasta(fasta_fp, smalla_fp)
        seq_cnt += 1


def write_header(fp, desc):
  """Given a freshly opened file, write the file header from the information given in the description file."""
  fp.seek(0)
  values = (__version__, desc['header']['species'].encode('utf8'), desc['header']['chromosome count'], desc['header']['ploidy of data'])
  fp.write(struct.pack(header_fmt, *values))


def append_index(fp, seq_info):
  fp.write(struct.pack(index_fmt, *(seq_info['chromosome number'], seq_info['chromosome copy'],
                                    seq_info['chromosome description'].encode('ascii'),
                                    seq_info['sequence id'],
                                    seq_info['start byte of sequence data'],
                                    seq_info['length of sequence'])))


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
  elif args['test']:
    import doctest
    doctest.testmod()
    exit(0)
  elif args['tofasta']:
    exit(0)
  else:
    with open(args['--descfile'], 'r') as f:
      desc = json.load(f)
    with open(args['--smalla'], 'w') as smalla_fp:
      create_smalla(desc, args['--data_dir'], smalla_fp)
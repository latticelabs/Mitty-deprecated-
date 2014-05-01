"""Given a ref_seq (smalla format) and a VCF file (indexed with tabix) construct the mutated genome as a smalla file
along with the POS index arrays. See Readme.md for detailed descriptions of the algorithms used by
vcf2seq.py and reads.py to generate reads with correct pos and CIGAR encoded (which can then be used by cheata.py to
generate perfect alignments and can be used to check alignment performance of other tools).

Usage:
vcf2seq <ref_seq>  <var_seq> <chrom> <vcf_file> [--block_size=BS] [-v]

Options:
  ref_seq            reference sequence in smalla format
  var_seq            output (variant) sequence (will be written in smalla format)
                     In addition .heada, .pos and .diffpos files are written
  chrom              Chromosome (Needed when pulling variants from indexed VCF file)
  vcf_file           indexed VCF file
  --block_size=BS    How many bases to process at one time [default: 1000000]
  -v                 Verbose (print debugging messages)

Notes:
1. If the VCF has alternate alleles only the first one is used when generating the variant sequence
2. Run tests by calling as vcf2seq test or vcf2seq test -v

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
__version__ = '0.2.0'
import docopt
import mmap  # For random, disk based access of reference .smalla files
import struct  # For writing pos and dpos data in binary format (need struct.pack)
import vcf
import logging

logger = logging.getLogger(__name__)

def block_copy_seq(ref_seq, start, stop, f_vseq, f_vseq_pos, block_size=1000):
  """This is a straight copy. It starts and ends at conserved sections.

  Try a complete copy
  >>> import io, struct; \
  ref_seq = 'ACTGACTG'; \
  f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  block_copy_seq(ref_seq, 0, 8, f_vseq, f_vseq_pos)
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  ACTGACTG
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('8I',f_vseq_pos.read(4*8))
  (1, 2, 3, 4, 5, 6, 7, 8)

  Copy just a section - note that the pos array should still contain the correct
  >>> f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  block_copy_seq(ref_seq, 3, 6, f_vseq, f_vseq_pos)
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  GAC
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('3I',f_vseq_pos.read(4*3))
  (4, 5, 6)


  Copy just a section with smallest block size possible
  >>> f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  block_copy_seq(ref_seq, 3, 6, f_vseq, f_vseq_pos, block_size=1)
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  GAC
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('3I',f_vseq_pos.read(4*3))
  (4, 5, 6)


  """
  this_start = start
  while this_start < stop:
    this_stop = min(this_start + block_size, stop)
    f_vseq.write(ref_seq[this_start:this_stop])
    f_vseq_pos.write(struct.pack('{:d}I'.format(stop-start), *[n+1 for n in range(start, stop)]))  # +1 because of
                                                                                                 # 1-indexing
    logger.debug('{:d}% done'.format(int(100.0 * this_stop / len(ref_seq))))
    this_start = this_stop


def handle_variant(variant, f_vseq, f_vseq_pos):
  """Write the variant to the mutated sequence and fill out the pos array

  Test with SNP
  >>> import io, vcf, struct; \
  ref_seq = 'ACTGACTG'; \
  pos = 2; \
  ref = 'C'; \
  alt = [vcf.model._Substitution('T')]; \
  variant = vcf.model._Record('1', pos, '.', ref, alt, 100, None, None, None, None, None); \
  f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos)
  2
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  T
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('I',f_vseq_pos.read(4))
  (2,)

  Test with delete
  >>> ref_seq = 'ACTGACTG'; \
  pos = 2; \
  ref = 'CTG'; \
  alt = [vcf.model._Substitution('C')]; \
  variant = vcf.model._Record('1', pos, '.', ref, alt, 100, None, None, None, None, None); \
  f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos)
  4
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  C
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('I',f_vseq_pos.read(4))
  (2,)


  Test with delete and .
  >>> ref_seq = 'ACTGACTG'; \
  pos = 8; \
  ref = 'G'; \
  alt = [vcf.model._Substitution('.')]; \
  variant = vcf.model._Record('1', pos, '.', ref, alt, 100, None, None, None, None, None); \
  f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos)
  8
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  <BLANKLINE>
  >>> _ = f_vseq_pos.seek(0,2); print f_vseq.tell()
  0


  Test with delete and .
  >>> ref_seq = 'ACTGACTG'; \
  pos = 4; \
  ref = 'GACTG'; \
  alt = [vcf.model._Substitution('.')]; \
  variant = vcf.model._Record('1', pos, '.', ref, alt, 100, None, None, None, None, None); \
  f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos)
  8
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  <BLANKLINE>
  >>> _ = f_vseq_pos.seek(0,2); print f_vseq.tell()
  0


  Test with insert
  >>> ref_seq = 'ACTGACTG'; \
  pos = 2; \
  ref = 'C'; \
  alt = [vcf.model._Substitution('CGGG')]; \
  variant = vcf.model._Record('1', pos, '.', ref, alt, 100, None, None, None, None, None); \
  f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos)
  2
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  CGGG
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('4I',f_vseq_pos.read(4*4))
  (2, 3, 3, 3)


  Test with insert
  >>> ref_seq = 'ACTGACTG'; \
  pos = 8; \
  ref = 'G'; \
  alt = [vcf.model._Substitution('GTTT')]; \
  variant = vcf.model._Record('1', pos, '.', ref, alt, 100, None, None, None, None, None); \
  f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos)
  8
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  GTTT
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('4I',f_vseq_pos.read(4*4))
  (8, 9, 9, 9)

  See Readme.md for details of algorithm
  """
  alt = variant.ALT[0].sequence
  if alt == '.': alt = ''
  ref = variant.REF
  if ref == '.': ref = ''
  f_vseq.write(alt)    # Copy over ALT
  new_ref_pos = len(ref) + variant.POS - 1  # Advance along ref_seq
                                            # -1 because POS is 1-indexed, we are 0-indexed internally
  if len(alt) > 0:
    if len(ref) > 0:
      f_vseq_pos.write(struct.pack('I', variant.POS))  # The original base
    if len(alt) > 1:
      f_vseq_pos.write(struct.pack('{:d}I'.format(len(alt) - len(ref)), *[new_ref_pos + 1] * (len(alt) - len(ref))))
  return new_ref_pos


def assemble_sequence(ref_seq, variants, f_vseq, f_vseq_pos, block_size=1000):
  """Given a ref_seq and list of variants, generate the var_seq and pos arrays.
  Inputs:
    ref_seq       - the reference sequence,
    variants      - any object that has a .next operator that will yield model._Record. This can be a list, or a generator
                    like vcf.Reader()
    f_vseq        - file like device for storing the mutated seq
    f_vseq_pos    - file like device for storing pos
    block_size    - how many bases we copy at a time. Only comes into play if the gap between variants is larger than
                    this

  Test with SNP, insert and delete
  >>> import io, vcf, struct; \
  variants = []; \
  ref_seq = 'ACTGACTG'; \
  pos = 2; \
  ref = 'C'; \
  alt = [vcf.model._Substitution('T')]; \
  variants.append(vcf.model._Record('1', pos, '.', ref, alt, 100, None, None, None, None, None)); \
  pos = 4; \
  ref = 'G'; \
  alt = [vcf.model._Substitution('GT')]; \
  variants.append(vcf.model._Record('1', pos, '.', ref, alt, 100, None, None, None, None, None)); \
  pos = 6; \
  ref = 'CTG'; \
  alt = [vcf.model._Substitution('.')]; \
  variants.append(vcf.model._Record('1', pos, '.', ref, alt, 100, None, None, None, None, None)); \
  f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  assemble_sequence(ref_seq, variants, f_vseq, f_vseq_pos);
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  ATTGTA
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('6I',f_vseq_pos.read(4*6))
  (1, 2, 3, 4, 5, 5)

  Now test the same thing, but with the smallest blocks possible. The answer should not change
  >>> f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  assemble_sequence(ref_seq, variants, f_vseq, f_vseq_pos, block_size=1);
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  ATTGTA
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('6I',f_vseq_pos.read(4*6))
  (1, 2, 3, 4, 5, 5)
  """
  copy_start = 0  # Start at the beginning
  for variant in variants:
    # Copy up to this variant
    block_copy_seq(ref_seq, copy_start, variant.POS - 1, f_vseq, f_vseq_pos, block_size)
    # -1 because we are zero indexed
    copy_start = handle_variant(variant, f_vseq, f_vseq_pos)

  # Now copy over any residual
  block_copy_seq(ref_seq, copy_start, len(ref_seq), f_vseq, f_vseq_pos, block_size)


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

  chrom = args['<chrom>']
  vcf_reader = vcf.Reader(filename=args['<vcf_file>'])
  ref_seq_fname = args['<ref_seq>']
  var_seq_fname = args['<var_seq>']
  var_seq_pos_fname = args['<var_seq>'] + '.pos'

  with open(ref_seq_fname, 'r+b') as f_rseq, open(var_seq_fname, 'w') as f_vseq,\
       open(var_seq_pos_fname, 'w') as f_vseq_pos:
    with open(var_seq_fname + '.heada', 'w') as fheada:  # Copy the header over
      fheada.write(open(ref_seq_fname + '.heada', 'r').readline() + ' (Mutated)')

    ref_seq = mmap.mmap(f_rseq.fileno(), 0)
    block_size = int(args['--block_size'])  # This is how many bases we copy at a time. Only comes into play if the
                                            # gap between variants is larger than this
    assemble_sequence(ref_seq, vcf_reader, f_vseq, f_vseq_pos, block_size)
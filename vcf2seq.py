"""Given a ref_seq (smalla format) and a VCF file (indexed with tabix) construct the mutated genome as a smalla file
along with the POS index arrays. See Readme.md for detailed descriptions of the algorithms used by
vcf2seq.py and reads.py to generate reads with correct pos and CIGAR encoded (which can then be used by cheata.py to
generate perfect alignments and can be used to check alignment performance of other tools).

Usage:
vcf2seq <ref_seq>  <var_seq> <chrom> <vcf_file> [--ploidy=PL] [--block_size=BS] [-v]

Options:
  ref_seq            reference sequence in smalla format
  var_seq            output (variant) sequence (will be written in smalla format)
                     In addition .heada, .pos and .diffpos files are written
  chrom              Chromosome (Needed when pulling variants from indexed VCF file)
  vcf_file           indexed VCF file
  --ploidy=PL        Use genotype information in VCF file to create polyploid sequences
                     If this is 1 (default) we ignore any genotype information in the VCF file and only produce one
                     .smalla mutated sequence. If this is set to 2 (or 3 or so on) we generate that many sequences
                     (numbered _1 _2 and so on) representing each strand of the sequence [default: 1].
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


def handle_variant(variant, f_vseq, f_vseq_pos, strand_no=None):
  """Write the variant to the mutated sequence and fill out the pos array

  Some setup
  >>> import io, vcf, struct
  >>> def add_GT(v, gt):
  ...   md=vcf.model.make_calldata_tuple('GT')(gt)
  ...   call=vcf.model._Call(v,'sample',md)
  ...   v.samples = [call]
  ...   v._sample_indexes = {'sample': 0}

  Test with SNP: ignore zygosity
  >>> ref_seq = 'ACTGACTG'; \
  pos = 2; \
  ref = 'C'; \
  alt = [vcf.model._Substitution('T')]; \
  variant = vcf.model._Record('1', pos, '.', ref, alt, 100, None, None, None, None, None); \
  add_GT(variant, '1/1'); \
  f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos)
  2
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  T
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('I',f_vseq_pos.read(4))
  (2,)

  Test with SNP: homo, strand 0
  >>> f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos, 0)
  2
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  T
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('I',f_vseq_pos.read(4))
  (2,)


  Test with SNP: homo, strand 1
  >>> f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos, 1)
  2
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  T
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('I',f_vseq_pos.read(4))
  (2,)



  Test with SNP: het, strand 0
  >>> add_GT(variant, '0/1'); \
  f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos, 0)
  1
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  <BLANKLINE>
  >>> _ = f_vseq_pos.seek(0,2); print f_vseq.tell()
  0

  Test with SNP: het, strand 1
  >>> f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos, 1)
  2
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  T
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('I',f_vseq_pos.read(4))
  (2,)


  Test with delete: ignore zygosity
  >>> ref_seq = 'ACTGACTG'; \
  pos = 2; \
  ref = 'CTG'; \
  alt = [vcf.model._Substitution('C')]; \
  variant = vcf.model._Record('1', pos, '.', ref, alt, 100, None, None, None, None, None); \
  add_GT(variant, '1/0'); \
  f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos)
  4
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  C
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('I',f_vseq_pos.read(4))
  (2,)

  Test with same delete, strand 1 (same as REF)
  >>> f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  print handle_variant(variant, f_vseq, f_vseq_pos, 1)
  1
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  <BLANKLINE>
  >>> _ = f_vseq_pos.seek(0,2); print f_vseq.tell()
  0


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
  var_type = '1' if strand_no is None else variant.genotype('sample').data.GT.split('/')[strand_no]
  # 0 means REF 1 means ALT. If we don't specify a strand number it means we don't care about Zygosity and will always
  # take the ALT
  if var_type == '1':
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
  else:
    new_ref_pos = variant.POS - 1  # Keep us here - we didn't implement this variant and we should keep copying
                                   # -1 because POS is 1-indexed, we are 0-indexed internally
  return new_ref_pos


def assemble_sequence(ref_seq, variants, f_vseq, f_vseq_pos, strand_no=None, block_size=1000):
  """Given a ref_seq and list of variants, generate the var_seq and pos arrays.
  Inputs:
    ref_seq       - the reference sequence,
    variants      - any object that has a .next operator that will yield model._Record. This can be a list, or a
                    generator like vcf.Reader()
    f_vseq        - file like device for storing the mutated seq
    f_vseq_pos    - file like device for storing pos
    strand_no     - If None, we ignore the GT information in the VCF records. If 0,1,2 ... we pick the value of the
                    relevant GT string to determine if the variant is to be applied or not.
    block_size    - how many bases we copy at a time. Only comes into play if the gap between variants is larger than
                    this

  Note: The way heterozygous genotyping works is that we put a number into strand_no and we look to the GT (Genotype)
  information in the VCF entry. Depending on the value of the relevant genotype (currently 0 or 1 only) we either write
  out the variant or we skip it.

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
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('7I',f_vseq_pos.read(4*7))
  (1, 2, 3, 4, 5, 5, 9)

  Now test the same thing, but with the smallest blocks possible. The answer should not change
  >>> f_vseq = io.BytesIO(); \
  f_vseq_pos = io.BytesIO(); \
  assemble_sequence(ref_seq, variants, f_vseq, f_vseq_pos, block_size=1);
  >>> _ = f_vseq.seek(0); print f_vseq.read()
  ATTGTA
  >>> _ = f_vseq_pos.seek(0); print struct.unpack('7I',f_vseq_pos.read(4*7))
  (1, 2, 3, 4, 5, 5, 9)
  """
  copy_start = 0  # Start at the beginning
  for variant in variants:
    # Copy up to this variant
    block_copy_seq(ref_seq, copy_start, variant.POS - 1, f_vseq, f_vseq_pos, block_size)
    # -1 because we are zero indexed
    copy_start = handle_variant(variant, f_vseq, f_vseq_pos, strand_no)

  # Now copy over any residual
  block_copy_seq(ref_seq, copy_start, len(ref_seq), f_vseq, f_vseq_pos, block_size)

  # Don't forget to add the 'tail': simply an additional, imaginary, base for housekeeping
  f_vseq_pos.write(struct.pack('I', len(ref_seq) + 1))


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

  chrom = args['<chrom>']  # TODO: use this info in reading the VCF
  ploidy = int(args['--ploidy'])
  ref_seq_fname = args['<ref_seq>']
  var_seq_fname_prefix = args['<var_seq>']
  #var_seq_pos_fname = args['<var_seq>'] + '.pos'
  with open(ref_seq_fname, 'rb') as f_rseq:
    f_vseq_l = [open('{:s}_{:d}.smalla'.format(var_seq_fname_prefix, n), 'w') for n in range(ploidy)]
    f_vseq_pos_l = [open('{:s}_{:d}.smalla.pos'.format(var_seq_fname_prefix, n), 'w') for n in range(ploidy)]
    header = open(ref_seq_fname + '.heada', 'r').readline() + ' (Mutated)'
    for n in range(ploidy): open('{:s}_{:d}.smalla.heada'.format(var_seq_fname_prefix, n), 'w').write(header)
    ref_seq = mmap.mmap(f_rseq.fileno(), 0, access=mmap.ACCESS_READ)
    block_size = int(args['--block_size'])  # This is how many bases we copy at a time. Only comes into play if the
                                            # gap between variants is larger than this
    for n, (f_vseq, f_vseq_pos) in enumerate(zip(f_vseq_l, f_vseq_pos_l)):
      vcf_reader = vcf.Reader(filename=args['<vcf_file>'])
      logger.debug('Writing strand {:d}'.format(n))
      strand_no = None if ploidy == 1 else n  # If ploidy is 1 we ignore the genotype info in the VCF if any
      assemble_sequence(ref_seq, vcf_reader, f_vseq, f_vseq_pos, strand_no=strand_no, block_size=block_size)
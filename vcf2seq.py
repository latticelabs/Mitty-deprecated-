"""Given a ref_seq (smalla format) and a VCF file (indexed with tabix) construct the mutated genome as a smalla file
along with the POS index arrays. See Readme.md for detailed descriptions of the algorithms used by
vcf2seq.py and reads.py to generate reads with correct pos and CIGAR encoded (which can then be used by cheata.py to
generate perfect alignments and can be used to check alignment performance of other tools).

Usage:
vcf2seq --ref=REF  --vcf=VCF  --var=SAMP  [-v]

Options:
  --ref=REF          reference genome sequence
  --vcf=VCF          vcf file
  --var=SAMP         output genome file (diploid)
  -v                 Verbose (print debugging messages)

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
__version__ = '0.3.0'
import h5py
import docopt
import numpy
import struct  # For writing pos and dpos data in binary format (need struct.pack)
import vcf
import logging

logger = logging.getLogger(__name__)
pos_null = numpy.empty((0,), dtype='u4')


def handle_variant(variant):
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
  #alt = variant.ALT[0].sequence
  #if alt == '.': alt = ''
  #ref = variant.REF
  #if ref == '.': ref = ''
  alt = variant.ALT[0].sequence if variant.ALT[0] is not None else ''
  ref = variant.REF or ''
  ptr_adv = len(ref) + 1  # You need to be one beyond ref
  pos_alt = numpy.arange(variant.POS, variant.POS + len(ref), dtype='u4')  # We might have ref bases in alt
  if len(alt) > len(ref):  # This was an insertion
    pos_alt = numpy.concatenate((pos_alt, numpy.ones(len(alt) - len(ref), dtype='u4') * (variant.POS + len(ref))))
  if variant.heterozygosity > 0:  # This is heterozygous
    if variant.samples[0].gt_nums[0] == '1':
      return [alt, ''], [pos_alt, pos_null], [ptr_adv, 0]
    else:
      return ['', alt], [pos_null, pos_alt], [0, ptr_adv]
  else:
    return [alt, alt], [pos_alt, pos_alt], [ptr_adv, ptr_adv]


def assemble_sequences(ref_seq, reader):
  copy = ['', '']
  pos = [pos_null, pos_null]
  pointer = [0, 0]  # Which part of the reference are we reading, for each copy of the chromosome
  for variant in reader:
    for n in [0, 1]:
      # Copy up to this variant
      copy[n] += ref_seq[pointer[n]:variant.POS - 1]
      pos[n] = numpy.concatenate((pos[n], numpy.arange(pointer[n] + 1, variant.POS, dtype='u4')))
      pointer[n] = variant.POS - 1

    # Handle variants
    var_seq, var_pos, ptr_adv = handle_variant(variant)
    for n in [0, 1]:
      copy[n] += var_seq[n]
      pos[n] = numpy.concatenate((pos[n], var_pos[n]))
      pointer[n] += ptr_adv[n]

  for n in [0, 1]:
    # Now copy over any residual
    copy[n] += ref_seq[pointer[n]:]
    pos[n] = numpy.concatenate((pos[n], numpy.arange(pointer[n] + 1, len(ref_seq) + 1, dtype='u4')))
    # Don't forget to add the 'tail': simply an additional, imaginary, base for housekeeping

  return copy, pos


def main():
  vcf_reader = vcf.Reader(filename=args['--vcf'])
  with h5py.File(args['--ref'], 'r') as ref_fp, h5py.File(args['--var'], 'w') as var_fp:
    for chrom in ref_fp['sequence']:
      ref_seq = ref_fp['sequence/{:d}/1'.format(chrom)][:].tostring()  # Very cheap operation
      copy, pos = assemble_sequences(ref_seq, vcf_reader.fetch(chrom=chrom, start=0, end=len(ref_seq)))
      for n, (this_copy, this_pos) in enumerate(zip(copy, pos)):
        var_fp.create_dataset('sequence/{:d}/{:d}'.format(chrom, n+1), data=this_copy)
        var_fp.create_dataset('pos/{:d}/{:d}'.format(chrom, n+1), data=this_pos)


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  main(cmd_args)
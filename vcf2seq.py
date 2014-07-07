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
import vcf
import logging

logger = logging.getLogger(__name__)
pos_null = numpy.empty((0,), dtype='u4')


def handle_variant(variant):
  """Write the variant to the mutated sequence and fill out the pos array"""
  alt = variant.ALT[0].sequence if variant.ALT[0] is not None else ''
  ref = variant.REF or ''
  ptr_adv = len(ref)  # When we get to this function, our pointer is sitting at POS
  pos_alt = numpy.arange(variant.POS, variant.POS + min(len(alt), len(ref)), dtype='u4')  # We might have ref bases in alt
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
      pointer[n] += max(0, variant.POS - 1 - pointer[n])

    # Handle variants
    var_seq, var_pos, ptr_adv = handle_variant(variant)
    for n in [0, 1]:
      copy[n] += var_seq[n]
      pos[n] = numpy.concatenate((pos[n], var_pos[n]))
      pointer[n] += ptr_adv[n]

  for n in [0, 1]:
    # Now copy over any residual
    copy[n] += ref_seq[pointer[n]:]
    pos[n] = numpy.concatenate((pos[n], numpy.arange(pointer[n] + 1, len(ref_seq) + 2, dtype='u4')))
    # Don't forget to add the 'tail': simply an additional, imaginary, base for housekeeping

  return copy, pos


def main(args):
  vcf_reader = vcf.Reader(filename=args['--vcf'])
  with h5py.File(args['--ref'], 'r') as ref_fp, h5py.File(args['--var'], 'w') as var_fp:
    for chrom in [int(c) for c in ref_fp['sequence']]:  #h5 keys are unicode
      # We assume the reference is haploid
      ref_seq = ref_fp['sequence/{:d}/1'.format(chrom)][:].tostring()  # Very cheap operation
      try:
        copy, pos = assemble_sequences(ref_seq, vcf_reader.fetch(chrom=chrom, start=0, end=len(ref_seq)))
      except KeyError:  # No mutations for this chromosome
        copy = [ref_seq, ref_seq]  # Same as reference
        pos = [None, None]
      for n, (this_copy, this_pos) in enumerate(zip(copy, pos)):
        var_fp.create_dataset('sequence/{:d}/{:d}'.format(chrom, n+1), data=numpy.fromstring(this_copy, dtype='u1'))
        if this_pos is not None:
          var_fp.create_dataset('pos/{:d}/{:d}'.format(chrom, n+1), data=this_pos)


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  main(cmd_args)
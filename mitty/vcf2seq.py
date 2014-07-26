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
pos_null = numpy.empty((0,), dtype='u4')  # Convenient, used in handle_variant

var_type = {'snp': 1, 'indel': 2, 'sv': 3}  # Needed to code the variant type


def handle_variant(variant):
  """Write the variant to the mutated sequence and fill out the pos array"""
  alt = variant.ALT[0].sequence if variant.ALT[0] is not None else ''
  ref = variant.REF or ''
  # The vcf module converts a '.' into a Python None

  ptr_adv = len(ref)  # When we get to this function, our pointer is sitting at POS
  pos_alt = numpy.arange(variant.POS, variant.POS + min(len(alt), len(ref)), dtype='u4')  # We might have ref bases in alt
  if len(alt) > len(ref):  # This was an insertion
    pos_alt = numpy.concatenate((pos_alt, numpy.ones(len(alt) - len(ref), dtype='u4') * (variant.POS + len(ref))))
  if variant.heterozygosity > 0:  # This is heterozygous
    if variant.samples[0].gt_nums[0] == '1':
      return [alt, ''], [pos_alt, pos_null], [ptr_adv, 0], 1  # Het copy 1
    else:
      return ['', alt], [pos_null, pos_alt], [0, ptr_adv], 2  # Het copy 2
  else:
    return [alt, alt], [pos_alt, pos_alt], [ptr_adv, ptr_adv], 0  # Homozygous


def assemble_sequences(ref_seq, reader):
  copy = [[], []]
  pos = [[pos_null], [pos_null]]
  pointer = [0, 0]  # Which part of the reference are we reading, for each copy of the chromosome
  mutated_pointer = [0, 0]  # Where are we on the mutated sequence (needed for marking variants)
  variant_coordinates = []  # n x 2 x 2 (variants, chromosome copies, [start, stop])
  variant_codes = []  # n x 2  (variants, [code, het])

  for variant in reader:
    for n in [0, 1]:
      # Copy up to this variant
      copy[n].append(ref_seq[pointer[n]:variant.POS - 1])
      pos[n].append(numpy.arange(pointer[n] + 1, variant.POS, dtype='u4'))
      pointer[n] += max(0, variant.POS - 1 - pointer[n])
      mutated_pointer[n] += len(copy[n][-1])

    # Handle variants
    var_seq, var_pos, ptr_adv, het = handle_variant(variant)
    variant_codes.append([var_type[variant.var_type], het])
    this_vc = []
    for n in [0, 1]:
      copy[n].append(var_seq[n])
      pos[n].append(var_pos[n])
      pointer[n] += ptr_adv[n]

      this_vc.append([mutated_pointer[n], mutated_pointer[n] + len(copy[n][-1])])

      mutated_pointer[n] += len(copy[n][-1])

    variant_coordinates.append(this_vc)

    logger.debug('({:d}%, {:d}%)'.format(100*pointer[0]/len(ref_seq), 100*pointer[1]/len(ref_seq)))

  for n in [0, 1]:
    # Now copy over any residual
    copy[n].append(ref_seq[pointer[n]:])
    pos[n].append(numpy.arange(pointer[n] + 1, len(ref_seq) + 2, dtype='u4'))
    # Don't forget to add the 'tail': simply an additional, imaginary, base for housekeeping

  return (''.join(copy[0]), ''.join(copy[1])), (numpy.concatenate(pos[0]), numpy.concatenate(pos[1])), variant_coordinates, variant_codes


def main(args):
  vcf_reader = vcf.Reader(filename=args['--vcf'])
  with h5py.File(args['--ref'], 'r') as ref_fp, h5py.File(args['--var'], 'w') as var_fp:
    logger.debug('Writing to {:s}'.format(var_fp.filename))
    var_fp.attrs['species'] = ref_fp.attrs['species']
    var_fp.create_group('sequence')
    for chrom in [int(c) for c in ref_fp['sequence']]:  # h5 keys are unicode
      logger.debug('Assembling chromosome {:d}'.format(chrom))
      # We assume the reference is haploid
      ref_seq = ref_fp['sequence/{:d}/1'.format(chrom)][:].tostring()  # Very cheap operation
      chrom_grp = var_fp['sequence'].create_group(str(chrom))
      chrom_grp.attrs['seq_id'] = ref_fp['sequence/{:d}'.format(chrom)].attrs['seq_id']
      try:
        copy, pos, var_coords, var_codes = assemble_sequences(ref_seq, vcf_reader.fetch(chrom=chrom, start=0, end=len(ref_seq)))
      except KeyError:  # No mutations for this chromosome
        logger.debug('Chromosome {:d} is unaltered, copying reference'.format(chrom))
        copy = [ref_seq, ref_seq]  # Same as reference
        pos = [None, None]
        var_coords, var_codes = [None, None]
      for n, (this_copy, this_pos) in enumerate(zip(copy, pos)):
        dset = chrom_grp.create_dataset('{:d}'.format(n+1), data=numpy.fromstring(this_copy, dtype='u1'))
        dset.attrs['reference'] = True
        if this_pos is not None:  # We have variants
          dset.attrs['reference'] = False
          var_fp.create_dataset('pos/{:d}/{:d}'.format(chrom, n+1), data=this_pos)
      if var_coords:
        var_fp.create_dataset('variants/pos/{:d}'.format(chrom), data=numpy.array(var_coords, dtype='u4'))
        var_fp.create_dataset('variants/codes/{:d}'.format(chrom), data=numpy.array(var_codes, dtype='u1'))


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  main(cmd_args)
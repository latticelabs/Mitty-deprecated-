"""Given a ref_seq (smalla format) and a VCF file (indexed with tabix) construct the mutated genome as a smalla file

Usage:
vcf2seq <ref_seq>  <var_seq> <chrom> <vcf_file> [--block_size=BS] [-v]

Options:
  ref_seq            reference sequence in smalla format
  var_seq            output (variant) sequence
  chrom              Chromosome
  vcf_file           index VCF file
  --block_size=BS    Block size [default: 1000000]
  -v                 Verbose (print debugging messages)

Currently only handles SNPs

TODO: Write this properly

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
__version__ = '0.1.0'
import docopt
import mmap
import vcf
import logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  elif docopt.sys.argv[1] == 'test':
    import doctest

    doctest.testmod()
    exit()
  else:
    args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)
  logging.debug(args)

  vcf_reader = vcf.Reader(filename=args['<vcf_file>'])

  snp_count = 0

  #Load the ref-seq smalla file
  f = open(args['<ref_seq>'], 'r+b')
  ref_seq = mmap.mmap(f.fileno(), 0)
  ref_seq_len = len(ref_seq)
  block_size = int(args['--block_size'])
  with open(args['<var_seq>'], 'w') as fout:
    block_start = 0
    while block_start < ref_seq_len:
      logger.debug('{:d}% done'.format(int(100 * block_start / float(ref_seq_len))))
      block_stop = block_start + block_size
      this_ref_seq_block = ref_seq[block_start:block_stop]
      variants = {rec.POS: rec for rec in
                  vcf_reader.fetch(args['<chrom>'], block_start + 1, min(ref_seq_len, block_stop))}

      #Slow code only for SNPS, FIXME
      for n, s in enumerate(this_ref_seq_block):
        if n + block_start + 1 in variants:
          fout.write(variants[n + block_start + 1].ALT[0].sequence)  #TODO warning about handling only 1st variant
          snp_count += 1
          #logger.debug('SNP!')
        else:
          fout.write(s)
      block_start += block_size

  logger.debug('{:d} SNPs'.format(snp_count))
  f.close()
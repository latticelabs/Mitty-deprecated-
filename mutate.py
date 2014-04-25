"""This module contains functions that generate variants with reference to a reference genome. The module can be called
as a script as well. This is useful for creating test data for MGR algorithms/data formats. The script will output
VCF file(s).

Usage:
mutate [--chrom=CHR]  --ref=REF  [--vcf=VCF]  --paramfile=PFILE [--block_size=BS] [-v]

Options:
  --chrom=CHROM           The VCF file needs chromosome info.
                          This is also available to the simulation if it needs it [default: 1]
  --ref=REF               The reference sequence in smalla format
  --vcf=VCF               The output VCF file. If not specified, the vcf file goes to stdout
  --paramfile=PFILE       Name for parameter file
  --block_size=BS         Block size for operations. Adjust to match memory/resources of platform [default: 100000]
                          This governs how many variants are generated at a time before being dumped to disk.
  -v                      Dump detailed logger messages

Note: Running the code without any arguments will print this help string and exit

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""

__version__ = '0.2.0'

import sys
import os
import imp
import json
import mmap  # To memory map our smalla files
import docopt
import datetime
import pysam  # For the tabix functions
import logging

logger = logging.getLogger(__name__)


def write_vcf_header(file_handle, sim_date, argv, reference_filename):
  """Given a file handle, write out a suitable header to start the VCF file
  Inputs:
    file_handle       - an open file handle
    sim_date          - date in string format. Get's dumped into 'fileDate'
    argv              - get's dumped into the 'source' string of the vcf file
    reference_filename  - name of the reference, dumped into the 'reference' string

  Notes:
  1. 'filedate' is the date of the simulation
  2. 'source' contains the version of the mutate program used to generate the data
  3. 'reference' contains the name of the file entered as the reference file

  """
  file_handle.write(
    """##fileformat=VCFv4.1
##fileDate={:s}
##source=mutate.py {:s} ({:s})
##reference={:s}
#CHROM POS     ID        REF    ALT     QUAL FILTER INFO\n""".format(sim_date, __version__, argv, reference_filename)
  )


def write_vcf_mutations(file_handle, chrom, variants):
  """Given a mutator format dictionary write the mutations in VCF format into the file
  Inputs:
    file_handle   - handle of an opened text file. The output will be appended to this file.
    chrom         - chromosome number
    variants      - list of tuples (POS, REF, ALT) : standard format as returned by the variant plugins
  """
  for var in variants:
    # Need to add +1 because first base is 1 while we are using 0 indexing internally
    file_handle.write("{:s}\t{:d}\t.\t{:s}\t{:s}\t96\tPASS\t.\n".format(chrom, var[0] + 1, var[1], var[2]))


'''
Variant location and conflict resolution:

Say we have v variant types that we are trying to generate.

Independently generate v lists of candidate positions cp
Create an empty 'footprint' list

Start with the earliest cp and generate a corresponding variant. Does the variant's footprint conflict with a footprint
already in the 'footprint' list? If so, discard this variant.

Move on to the next closest variant position and continue


'''

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)
  logging.debug(args)

  params = json.load(open(args['--paramfile'], 'r')) #imp.load_source('params', args['--paramfile'], open(args['--paramfile'], 'r'))
  block_size = int(args['--block_size'])

  #Load the ref-seq smalla file
  fin = open(args['--ref'], 'r+b')
  ref_seq = mmap.mmap(fin.fileno(), 0)
  ref_seq_len = len(ref_seq)

  fout = sys.stdout if args['--vcf'] is None else open(args['--vcf'], 'w')

  plugin_dir = os.path.join(os.path.dirname(__file__), 'Plugins', 'Mutation')  # Thanks Nebojsa Tijanic!
  variant_generator = {}
  next_variant = {}
  misc = {}  # For now, holds a count of how many variants of each type
  for k in params.keys():
    model_fname = os.path.join(plugin_dir, params[k]['model'] + '_plugin.py')
    model = imp.load_source(k, model_fname, open(model_fname, 'r'))
    variant_generator[k] = model.variant(ref_seq=ref_seq, ref_seq_len=ref_seq_len, block_size=block_size, **params[k])
    next_variant[k] = next(variant_generator[k], None)  # Will return None if we have no more of these variants
    misc[k] = 0

  logger.debug('Input sequence has {:d} bases'.format(ref_seq_len))

  current_pos = 0  # Current position on sequence

  write_vcf_header(fout, datetime.datetime.now().isoformat(), docopt.sys.argv.__str__(), args['--ref'])

  #For now, ignoring footprint (only needed for distal variants)
  these_variants = []
  while current_pos < ref_seq_len:
    next_variant_pos = ref_seq_len
    next_variant_type = None
    # Find the earliest variant
    for k in params.keys():
      if next_variant[k] is not None:
        if next_variant[k][0] < next_variant_pos:
          next_variant_pos = next_variant[k][0]
          next_variant_type = k

    if next_variant_type is not None:  # We are still in business
      these_variants.append(next_variant[next_variant_type])  # Queue this variant
      misc[next_variant_type] += 1
      skip_to = next_variant[next_variant_type][3]  # The next variant has to come at or after this to avoid conflicts
      # Move all variant counters forward as needed
      for k in params.keys():
        while next_variant[k] is not None:
          if next_variant[k][0] < skip_to:
            next_variant[k] = next(variant_generator[k], None)
          else:
            break

    # Flush variants to disk
    if len(these_variants) == block_size:
      write_vcf_mutations(fout, args['--chrom'], these_variants)
      logger.debug('{:d}% done'.format(int(100.0 * current_pos / ref_seq_len)))
      these_variants = []

    current_pos = next_variant_pos

  # Flush any remaining variants to disk
  write_vcf_mutations(fout, args['--chrom'], these_variants)

  fin.close()
  fout.close()
  for k in params.keys():
    logger.debug('Generated {:d} {:s}s'.format(misc[k], k))

  logger.debug('Compressing and indexing VCF file')
  pysam.tabix_compress(args['--vcf'], args['--vcf'] + '.gz', force=True)
  pysam.tabix_index(args['--vcf'] + '.gz', force=True, preset='vcf')

  if args['--vcf'] is not None:
    with open(args['--vcf'] + '.info','w') as f:
      f.write('Command line\n-------------\n')
      f.write(json.dumps(args, indent=4))
      f.write('\n\nParameters\n------------\n')
      f.write(json.dumps(params, indent=4))

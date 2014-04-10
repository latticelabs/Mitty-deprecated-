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
import logging

logger = logging.getLogger(__name__)


# TODO: Use vcf tools for this?
def write_vcf_header(file_handle, sim_date, argv, reference_filename):
  """Given a file handle, write out a suitable header to start the VCF file
  Inputs:
    file_handle       - an open file handle

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


# TODO Update docs
def write_vcf_mutations(file_handle, chrom, variants):
  """Given a mutator format dictionary write the mutations in VCF format into the file
  Inputs:
    file_handle   - handle of an opened text file. The output will be appended to this file.
    mutations     - dictionary of lists in mutator format
  Notes:
  1. The mutator format is a dictionary of lists which makes it easy to drop the information into the VCF file.
     The dictionary keys correspond to the VCF columns. Each value is a list with however many elements we want the
     VCF file to have
  2. The ID is a madeup id that is guaranteed not to clash with any dbSNP id because it is a string unique to the
     mutate program
  """
  for var in variants:
    # Need to add +1 because first base is 1 while we are using 0 indexing internally
    file_handle.write("{:s}\t{:d}\t.\t{:s}\t{:s}\t96\tPASS\t.\n".format(chrom, var[0] + 1, var[1], var[2]))


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if args['-v'] else logging.WARNING
  logging.basicConfig(level=level)
  logging.debug(args)

  params = json.load(open(args['--paramfile'], 'r')) #imp.load_source('params', args['--paramfile'], open(args['--paramfile'], 'r'))

  #Load the ref-seq smalla file
  fin = open(args['--ref'], 'r+b')
  ref_seq = mmap.mmap(fin.fileno(), 0)
  ref_seq_len = len(ref_seq)

  fout = sys.stdout if args['--vcf'] is None else open(args['--vcf'], 'w')

  plugin_dir = os.path.join(os.path.dirname(__file__), 'Plugins', 'Mutation')  # Thanks Nebojsa Tijanic!
  model = {}
  for k in params.keys():
    model_fname = os.path.join(plugin_dir, params[k]['model'] + '_plugin.py')
    model[k] = imp.load_source(k, model_fname, open(model_fname, 'r'))

  prev_state = None
  logger.debug('Input sequence has {:d} bases'.format(ref_seq_len))
  block_size = int(args['--block_size'])
  block_start = 0
  variant_count_snp = 0
  write_vcf_header(fout, datetime.datetime.now().isoformat(), docopt.sys.argv.__str__(), args['--ref'])
  prev_state = {k: None for k in params.keys()}
  while block_start < ref_seq_len:
    this_ref_seq_block = ref_seq[block_start:block_start + block_size]
    these_variants = {}
    for k in params.keys():
      these_variants[k], prev_state[k] = \
        model[k].candidate_variants(chrom=args['--chrom'],
                                    ref_seq_len=ref_seq_len,
                                    ref_seq_block_start=block_start,
                                    ref_seq_block=this_ref_seq_block,
                                    prev_state=prev_state[k],
                                    **params[k])

      variant_count_snp += len(these_variants[k])  # FIXME

    #Need to resolve variants here

    for k in params.keys():
      write_vcf_mutations(fout, args['--chrom'], these_variants[k])

    block_start += block_size
    logger.debug('{:d}% done'.format(int(100.0 * min(block_start, ref_seq_len) / float(ref_seq_len))))

  fin.close()
  fout.close()
  logger.debug('Generated {:d} SNPs'.format(variant_count_snp))

  if args['--vcf'] is not None:
    with open(args['--vcf'] + '.info','w') as f:
      f.write('Command line\n-------------\n')
      f.write(json.dumps(args, indent=4))
      f.write('\n\nParameters\n------------\n')
      f.write(json.dumps(params, indent=4))

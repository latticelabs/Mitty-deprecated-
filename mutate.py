"""This module contains functions that generate variants with reference to a reference genome. The module can be called
as a script as well. This is useful for creating test data for MGR algorithms/data formats. The script will output
VCF file(s).

Usage:
mutate --paramfile=PFILE [--block_size=BS] [-v]

Options:
  --paramfile=PFILE       Name for parameter file
  --block_size=BS         Block size for operations. Adjust to match memory/resources of platform [default: 100000]
  -v                      Dump detailed logger messages

Note: Running the code without any arguments will print this help string and exit

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""

__version__ = '0.2.0'

import sys
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
  fin = open(params['ref'], 'r+b')
  ref_seq = mmap.mmap(fin.fileno(), 0)
  ref_seq_len = len(ref_seq)

  fout = sys.stdout if params['vcf'] is u'' else open(params['vcf'], 'w')

  model_fname = 'Plugins/Mutation/' + params['models']['snp']['model'] + '_plugin.py'  # TODO: join paths properly
  snp_model = imp.load_source('snps', model_fname, open(model_fname, 'r'))

  prev_state = None
  logger.debug('Input sequence has {:d} bases'.format(ref_seq_len))
  block_size = int(args['--block_size'])
  block_start = run_start = params['models']['snp']['start']
  run_stop = params['models']['snp']['stop']
  if run_stop == -1: run_stop = ref_seq_len
  write_vcf_header(fout, datetime.datetime.now().isoformat(), docopt.sys.argv.__str__(), params['ref'])
  while block_start < run_stop:
    logger.debug('{:d}% done'.format(int(100 * (block_start - run_start) / float(run_stop - run_start))))
    this_ref_seq_block = ref_seq[block_start:block_start + block_size]
    snp_variants, prev_state = \
      snp_model.candidate_variants(chrom=params['chrom'],
                                   start_loc=block_start,
                                   ref_seq=this_ref_seq_block,
                                   prev_state=prev_state,
                                   **params['models']['snp']['args'])
    #Need to resolve variants here
    write_vcf_mutations(fout, params['chrom'], snp_variants)
    block_start += block_size

  fin.close()
  fout.close()
  sidecar_fname = params['vcf'] if len(params['vcf']) else 'temporary.vcf'  # Use this name when we pipe to stdout
  with open(sidecar_fname + '.info','w') as f:  # TODO: Use os.join
    f.write('Command line\n-------------\n')
    f.write(json.dumps(args, indent=4))
    f.write('\n\nParameters\n------------\n')
    f.write(json.dumps(params, indent=4))

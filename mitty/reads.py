#!python
__cmd__ = """This program generates reads from a genome (either a reference genome or a sample genome).
The read characteristics are governed by the chosen read plugin.

Commandline::

  Usage:
    reads  --pfile=PFILE  [--coverage_per_block=CPB] [-v|-V]
    reads list
    reads explain (parameters|model <model_name>)

  Options:
    --pfile=PFILE   Name for parameter file
    --coverage_per_block=CPB    Compute reads in blocks, each with this coverage [default: 1]
    list            List available models
    explain         Explain details about the indicated model
    parameters      Explain parameters
    model           Explain a model
    <model_name>    The model to explain
    -v              Dump detailed logger messages
    -V              Dump very detailed logger messages
"""

__param__ = """Parameter file example::

  {
    "files": {
      # An absolute path is left as is
      # a relative path is taken relative to the location of the *script*
      "reference_dir": "/Users/kghose/Data/hg38/",  # Use this if the reference consists of multiple .fa files in a directory
      "reference_file": "/Users/kghose/Data/hg38/hg38.fa.gz",  # Use this if reference is a single multi-fasta file
      "dbfile": "Out/test.db"  # Genomes database file. Leave out if taking reads from a reference, or from VCF file
      "input_vcf": "Out/test.vcf",  # Use this if using a VCF file. Leave out if taking reads from a reference, or from VCF file
      "output_prefix": "Out/reads"  # Output file name prefix
                                    # the reads will be called reads.fq and reads_c.fq if we call for corrupted reads too
    },
    "sample" : {
      "gen": 0,      # Use gen and serial for database
      "serial": 1,
      "name": "HN23456" # Use name for VCF file
    }
    "rng": {
      "master_seed": 1
    },
    "chromosomes": [1,2],               # List the chromosomes the reads should be taken from
    "variants_only": false,             # If true, reads will only come from the vicinity of variants
    "variant_window": 500,              # Only used if variants_only is true. Reads will be taken from this vicinity
    "corrupt": true,                    # If true, corrupted reads will also be written
    "coverage": 10,                     # Coverage
    "read_model": "simple_illumina",    # Model specific parameters, need to be under the key "model_params"
    "model_params": {
      "read_len": 100,          # length of each read
      "template_len_mean": 250, # mean length of template
      "template_len_sd": 30,    # sd of template length
      "max_p_error": 0.01,      # Maximum error rate at tip of read
      "k": 20                   # exponent
    }

  }
"""
__doc__ = __cmd__ + __param__
# We split this up so that when we print the help, we can print it as separate pages. It got cluttered when printed all
# together. We want this to be the docstring because that is a nice format for our generated documentation

import json
import os
import time
import sys

import numpy as np
import docopt

import mitty.lib
import mitty.lib.reads as lib_reads
import mitty.lib.io as mio
import mitty.lib.db as mdb
import mitty.lib.variants as vr

import logging
logger = logging.getLogger(__name__)


def cli():  # pragma: no cover
  """Serves as entry point for scripts"""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__cmd__, ['-h'])
  else:
    cmd_args = docopt.docopt(__cmd__)  # Version string?

  level = logging.DEBUG if cmd_args['-V'] else logging.WARNING
  logging.basicConfig(level=level)

  if cmd_args['-v']:
    logger.setLevel(logging.DEBUG)

  if cmd_args['list']:
    print_list(cmd_args)
  elif cmd_args['explain']:
    explain(cmd_args)
  else:
    run(cmd_args)


def print_list(cmd_args):
  print('\nAvailable models\n----------------')
  for name, mod_name in mitty.lib.discover_all_reads_plugins():
    print('- {:s} ({:s})\n'.format(name, mod_name))


def explain(cmd_args):
  if cmd_args['parameters']:
    print(__param__)
  elif cmd_args['model']:
    explain_read_model(cmd_args['<model_name>'])


# TODO update this when we use classes for read models
def explain_read_model(name):
  try:
    mod = mitty.lib.load_reads_plugin(name)
  except ImportError:
    print('No model named {:s}'.format(name))
    return
  try:
    print(mod._description)
  except AttributeError:
    print('No help for model "{:s}" available'.format(name))
  return


def run(cmd_args):
  """The main read generating loop

  :param cmd_args: parameters as parsed by doc_opt
  """
  base_dir = os.path.dirname(cmd_args['--pfile'])     # Other files will be with respect to this
  params = json.load(open(cmd_args['--pfile'], 'r'))

  null_reads = True
  if 'dbfile' in params['files']:
    pop_db_name = mitty.lib.rpath(base_dir, params['files']['dbfile'])
    conn = mdb.connect(db_name=pop_db_name)
    null_reads = False
  elif 'input_vcf' in params:
    raise(NotImplementedError, 'Reading from VCF not implemented yet')

  if not null_reads:
    gen = params['sample']['gen']
    serial = params['sample']['serial']

  ref = mio.Fasta(multi_fasta=mitty.lib.rpath(base_dir, params['files'].get('reference_file', None)),
                  multi_dir=mitty.lib.rpath(base_dir, params['files'].get('reference_dir', None)),
                  persistent=False)
  master_seed = int(params['rng']['master_seed'])
  assert 0 < master_seed < mitty.lib.SEED_MAX

  seed_rng = np.random.RandomState(seed=master_seed)

  chromosomes = params['chromosomes']
  #TODO: make sure cpb is sane
  coverage = float(params['coverage'])
  coverage_per_block = float(cmd_args['--coverage_per_block'])
  blocks_to_do = int(0.5 * coverage / coverage_per_block)  # Coverage reduced by 2 because we have two chromosome copies
  coverage_per_block = coverage / blocks_to_do
  read_model = mitty.lib.load_reads_plugin(params['read_model']).Model(**params['model_params'])
  corrupt = bool(params['corrupt'])

  t0 = time.time()
  first_serial_no = 0
  fname_prefix = mitty.lib.rpath(base_dir, params['files']['output_prefix'])
  fastq_fp = open(fname_prefix + '.fq', 'w')
  fastq_c_fp = open(fname_prefix + '_c.fq', 'w') if corrupt else None

  total_blocks = len(chromosomes) * 2 * blocks_to_do

  blocks_done = 0
  for ch in chromosomes:
    if not null_reads:
      ml = mdb.load_master_list(conn, ch)
      chrom = mdb.load_sample(conn, gen, serial, ch)
    else:
      ml, chrom = vr.VariantList(), []  # Need a dummy variant list for nulls
    for cpy in [0, 1]:
      seq, variant_waypoints = lib_reads.expand_sequence(ref[ch], ml, chrom, cpy)
      seq_c = mitty.lib.string.translate(seq, mitty.lib.DNA_complement)
      for blk in range(blocks_to_do):
        if not params['variants_only']:
          reads, paired = read_model.get_reads(seq, seq_c, coverage=coverage_per_block, corrupt=corrupt)
        else:
          raise(NotImplementedError, 'Reads from variants only is being developed')
        pos, cigars = lib_reads.roll_cigars(variant_waypoints, reads)
        first_serial_no = write_reads_to_file(fastq_fp, fastq_c_fp, reads, paired, pos, cigars, ch, cpy, first_serial_no)
        blocks_done += 1
        progress_bar('Generating reads ', f=float(blocks_done)/total_blocks, cols=80)
  print('')
  t1 = time.time()
  logger.debug('Took {:f}s to write {:d} templates'.format(t1 - t0, first_serial_no))


def reads_from_variants_only(seq, seq_c, ml, chrom, cpy, variant_window, coverage_per_block, corrupt, seed):
  """To do

  :param seq:      forward sequence
  :param seq_c:    complement sequence
  :param ml:       master list of variants
  :param chrom:    list of pointers to master list
  :param cpy:      which copy of chromosome we are considering
  :param variant_window: how many bases before and after variant should we include
  :param coverage_per_block: coverage level per batch of reads
  :param corrupt:  T/F generate corrupted reads too or not
  :param seed:     seed for the simulation
  :return:
  """
  pass


def write_reads_to_file(fastq_fp, fastq_c_fp, reads, paired, pos, cigars, ch, cc, first_serial_no):
  """
  :param fastq_fp:     file pointer to perfect reads file
  :param fastq_c_fp:   file pointer to corrupted reads file (None if no corrupted reads being written)
  :param reads:        reads recarray
  :param paired:       bool, are reads paired
  :param pos:          list of POS values for the reads
  :param cigars:       list of cigar values for the reads
  :param ch:           chromosome number
  :param cc:           chromosome copy
  :param first_serial_no: serial number of first template in this batch
  :return: next_serial_no: the serial number the next batch should start at
  """
  cntr = first_serial_no
  pr_seq, cr_seq, phred = reads['perfect_reads'], reads['corrupt_reads'], reads['phred']

  if paired:
    for n in xrange(0, reads.shape[0], 2):
      qname = 'r{:d}|{:d}:{:d}|{:d}:{:s}|{:d}:{:s}'.format(cntr, ch, cc, pos[n], cigars[n], pos[n + 1], cigars[n + 1])
      fastq_fp.write('@' + qname + '\n' + pr_seq[n] + '\n+\n' + '~' * len(pr_seq[n]) + '\n')
      fastq_fp.write('@' + qname + '\n' + pr_seq[n + 1] + '\n+\n' + '~' * len(pr_seq[n + 1]) + '\n')
      if fastq_c_fp is not None:
        fastq_c_fp.write('@' + qname + '\n' + cr_seq[n] + '\n+\n' + phred[n] + '\n')
        fastq_c_fp.write('@' + qname + '\n' + cr_seq[n + 1] + '\n+\n' + phred[n + 1] + '\n')
      cntr += 1
  else:
    for n in xrange(0, reads.shape[0]):
      qname = 'r{:d}|{:d}:{:s}'.format(cntr, pos[n], cigars[n])
      fastq_fp.write('@' + qname + '\n' + pr_seq[n] + '\n+\n' + '~' * len(pr_seq[n]) + '\n')
      if fastq_c_fp is not None:
        fastq_c_fp.write('@' + qname + '\n' + cr_seq[n] + '\n+\n' + phred[n] + '\n')
      cntr += 1

  return cntr


def progress_bar(title, f, cols):
  """Draw a nifty progress bar.
  '\r' trick from http://stackoverflow.com/questions/15685063/print-a-progress-bar-processing-in-python

  :param title: leading text to print
  :param f:     fraction completed
  :param cols:  how many columns wide should the bar be
  """
  x = int(f * cols + 0.5)
  sys.stdout.write('\r' + title + '[' + '.' * x + ' ' * (cols - x) + ']')
  sys.stdout.flush()


if __name__ == '__main__':
  cli()
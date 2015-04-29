#!python
__cmd__ = """This program generates reads from a genome (either a reference genome or a sample genome).
The read characteristics are governed by the chosen read plugin.

Commandline::

  Usage:
    reads generate <pfile>  [-v|-V] [-p]
    reads list
    reads explain (parameters|model (all|<model_name>))

  Options:
    <pfile>         Name for parameter file
    list            List available models
    explain         Explain details about the indicated model
    parameters      Explain parameters
    model           Explain a model
    all             Iterate over all models and explain them
    <model_name>    The model to explain
    -p              Show progress bar
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
      "output_prefix": "Out/reads", # Output file name prefix
                                    # the reads will be called reads.fq and reads_c.fq if we call for corrupted reads too
      "interleaved": true     # Set to false if you need separate files for each mate of a pair
                              # files will then be named reads_1.fq, reads_c_1.fq and reads_2.fq, reads_c_2.fq
                              # If the reads are actually not paired and you set interleaved to false you will get two
                              # files _1 and _2 and all the data will be in _1 only
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
    "coverage_per_block": 0.01          # For each block of the simulation we generate reads giving this coverage
                                        # In some simulations we will have too many reads to fit in memory and we must
                                        # write them out in blocks of smaller size.
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
    executor(cmd_args)


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
    print(mitty.lib.model_init_signature_string(mod.Model.__init__))
  except AttributeError:
    print('No help for model "{:s}" available'.format(name))
  return


def explain_all_read_models():
  for name, mod_name in mitty.lib.discover_all_reads_plugins():
    explain_read_model(name)


def generate_reads_loop(ref, conn=None, gen=None, serial=None, chromosomes=[], read_model=None,
                        variants_only=False, variant_window=None,
                        fastq_fp=[], fastq_c_fp=[],
                        actual_coverage_per_block=1.0, blocks_to_do=1,
                        master_seed=1,
                        progress_bar_func=None):
  """

  :param ref: reference sequence object
  :param conn: genome database connection. Leave None for null reads
  :param gen:    generation
  :param serial: serial number (leave both out for null reads)
  :param chromosomes: list of chromosomes to take reads from
  :param read_model: read model object
  :param variants_only: True if reads should be taken from regions neighboring variants only
  :param variant_window: Size of region (in bases) around variants for local reads
  :param fastq_fp:   List of two file pointers to perfect reads files. If interleaved, the pointers are same
  :param fastq_c_fp: List of two file pointers to corrupted reads file, leave None to skip computing corrupted reads.
                     If interleaved, the pointers are same
  :param actual_coverage_per_block: coverage of reads per computation block
  :param blocks_to_do: total blocks that need to be done per pass
  :param master_seed: rng seed
  :param progress_bar_func: function that draws a progress bar. Leave None to skip printing progress bar
  :return: total number of reads generated
  """
  first_serial_no = 0
  total_blocks = len(chromosomes) * 2 * blocks_to_do  # For the progress bar
  blocks_done = 0
  seed_rng = np.random.RandomState(seed=master_seed)
  corrupt = fastq_c_fp is not None

  if progress_bar_func: progress_bar_func('Generating reads ', f=float(blocks_done) / total_blocks, cols=80)
  for ch in chromosomes:
    if conn is not None:
      ml = mdb.load_master_list(conn, ch)
      chrom = mdb.load_sample(conn, gen, serial, ch)
    else:
      ml, chrom = vr.VariantList(), []  # Need a dummy variant list for nulls
    for cpy in [0, 1]:
      seq, variant_waypoints = lib_reads.expand_sequence(ref[ch]['seq'], ml, chrom, cpy)
      seq_c = mitty.lib.string.translate(seq, mitty.lib.DNA_complement)
      for blk in range(blocks_to_do):
        if not variants_only:
          reads, paired = read_model.get_reads(seq, seq_c, coverage=actual_coverage_per_block, corrupt=corrupt,
                                               seed=seed_rng.randint(0, mitty.lib.SEED_MAX))
        else:
          reads, paired = reads_from_variants_only(seq, seq_c, variant_waypoints, variant_window, read_model,
                                                   actual_coverage_per_block, corrupt, seed_rng)
        pos, cigars = lib_reads.roll_cigars(variant_waypoints, reads)
        first_serial_no = write_reads_to_file(fastq_fp[0], fastq_fp[1], fastq_c_fp[0], fastq_c_fp[1],
                                              reads, paired, pos, cigars, ch, cpy,
                                              first_serial_no)
        blocks_done += 1
        if progress_bar_func: progress_bar_func('Generating reads ', f=float(blocks_done) / total_blocks, cols=80)
  if progress_bar_func: print('')
  return first_serial_no


def executor(cmd_args):
  """Executor for read generation. Parses command line arguments, sets up file handles and then off we go!

  :param cmd_args: parameters as parsed by doc_opt
  """
  base_dir = os.path.dirname(cmd_args['<pfile>'])     # Other files will be with respect to this
  params = json.load(open(cmd_args['<pfile>'], 'r'))

  null_reads = True
  if 'dbfile' in params['files']:
    pop_db_name = mitty.lib.rpath(base_dir, params['files']['dbfile'])
    conn = mdb.connect(db_name=pop_db_name)
  elif 'input_vcf' in params:
    raise(NotImplementedError, 'Reading from VCF not implemented yet')
  else:
    conn = None
    logger.debug('Taking reads from reference')

  if conn is not None:
    gen = params['sample']['gen']
    serial = params['sample']['serial']
  else:
    gen, serial = None, None

  ref = mio.Fasta(multi_fasta=mitty.lib.rpath(base_dir, params['files'].get('reference_file', None)),
                  multi_dir=mitty.lib.rpath(base_dir, params['files'].get('reference_dir', None)),
                  persistent=False)
  master_seed = int(params['rng']['master_seed'])
  assert 0 < master_seed < mitty.lib.SEED_MAX

  chromosomes = params['chromosomes']

  coverage_per_chromosome = float(params['coverage']) / 2.0
  coverage_per_block = float(params['coverage_per_block'])

  blocks_to_do = int(coverage_per_chromosome / coverage_per_block)
  actual_coverage_per_block = coverage_per_chromosome / blocks_to_do  # Actual coverage has to be adjusted, because blocks_to_do is an int

  read_model = mitty.lib.load_reads_plugin(params['read_model']).Model(**params['model_params'])
  corrupt = bool(params['corrupt'])

  variants_only = params.get('variants_only', None)
  variant_window = int(params.get('variant_window', 200)) if variants_only else None

  t0 = time.time()
  fname_prefix = mitty.lib.rpath(base_dir, params['files']['output_prefix'])
  if params['files'].get('interleaved', True):
    fastq_fp = [open(fname_prefix + '.fq', 'w')] * 2
    fastq_c_fp = [open(fname_prefix + '_c.fq', 'w')] * 2 if corrupt else [None, None]
  else:  # Need two files separately
    fastq_fp = [open(fname_prefix + '_1.fq', 'w'), open(fname_prefix + '_2.fq', 'w')]
    fastq_c_fp = [open(fname_prefix + '_c_1.fq', 'w'), open(fname_prefix + '_c_2.fq', 'w')] if corrupt else [None, None]

  read_count = generate_reads_loop(ref=ref, conn=conn, gen=gen, serial=serial, chromosomes=chromosomes,
                                   read_model=read_model, variants_only=variants_only, variant_window=variant_window,
                                   fastq_fp=fastq_fp, fastq_c_fp=fastq_c_fp,
                                   actual_coverage_per_block=actual_coverage_per_block, blocks_to_do=blocks_to_do,
                                   master_seed=master_seed,
                                   progress_bar_func=mitty.lib.progress_bar if cmd_args['-p'] else None)
  t1 = time.time()
  logger.debug('Took {:f}s to write {:d} templates'.format(t1 - t0, read_count))


def reads_from_variants_only(seq, seq_c, variant_waypoints, variant_window, read_model, coverage, corrupt, seed_rng):
  """Wrapper around read function to get reads only from around variant locations

  :param seq:      forward sequence
  :param seq_c:    complement sequence
  :param variant_waypoints: as returned by expand_sequence
  :param variant_window: how many bases before and after variant should we include
  :param corrupt:  T/F generate corrupted reads too or not
  :param seed_rng:    rng for seed generation
  :return:
  """
  reads, paired = [], False
  for v in variant_waypoints:  # [1:-1]:
    start, stop = max(v[1] - variant_window, 0), min(v[1] + variant_window, len(seq))
    # v[1] is the pos of the variant in sequence coordinates (rather than ref coordinates)
    these_reads, paired = read_model.get_reads(seq, seq_c, start_base=start, end_base=stop, coverage=coverage, corrupt=corrupt, seed=seed_rng.randint(0, mitty.lib.SEED_MAX))
    reads += [these_reads]
  return np.concatenate(reads), paired


def write_reads_to_file(fastq_fp_1, fastq_fp_2, fastq_c_fp_1, fastq_c_fp_2,
                        reads, paired, pos, cigars, ch, cc, first_serial_no):
  """
  :param fastq_fp_1:   file pointer to perfect reads file
  :param fastq_fp_2:   file pointer to perfect reads file 2. Same as 1 if interleaving. None if not paired
  :param fastq_c_fp_1: file pointer to corrupted reads file. None if no corrupted reads being written
  :param fastq_c_fp_2: file pointer to corrupted reads file 2. Same as 1 if interleaving. None if no corrupted reads being written. None if not paired
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
  ro, pr_seq, cr_seq, phred = reads['read_order'], reads['perfect_reads'], reads['corrupt_reads'], reads['phred']

  if paired:
    for n in xrange(0, reads.shape[0], 2):
      qname = 'r{:d}|{:d}|{:d}|{:d}|{:d}|{:s}|{:d}|{:d}|{:s}'.format(cntr, ch, cc, ro[n], pos[n], cigars[n], ro[n + 1], pos[n + 1], cigars[n + 1])
      fastq_fp_1.write('@' + qname + ' /1\n' + pr_seq[n] + '\n+\n' + '~' * len(pr_seq[n]) + '\n')
      fastq_fp_2.write('@' + qname + ' /2\n' + pr_seq[n + 1] + '\n+\n' + '~' * len(pr_seq[n + 1]) + '\n')
      if fastq_c_fp_1 is not None:
        fastq_c_fp_1.write('@' + qname + ' /1\n' + cr_seq[n] + '\n+\n' + phred[n] + '\n')
        fastq_c_fp_2.write('@' + qname + ' /2\n' + cr_seq[n + 1] + '\n+\n' + phred[n + 1] + '\n')
      cntr += 1
  else:
    for n in xrange(0, reads.shape[0]):
      qname = 'r{:d}|{:d}|{:d}|{:d}|{:d}|{:s}'.format(cntr, ch, cc, ro[n], pos[n], cigars[n])
      fastq_fp_1.write('@' + qname + '\n' + pr_seq[n] + '\n+\n' + '~' * len(pr_seq[n]) + '\n')
      if fastq_c_fp_1 is not None:
        fastq_c_fp_1.write('@' + qname + '\n' + cr_seq[n] + '\n+\n' + phred[n] + '\n')
      cntr += 1

  return cntr


if __name__ == '__main__':
  cli()
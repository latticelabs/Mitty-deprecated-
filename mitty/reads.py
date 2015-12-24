#!python
__param__ = """Parameter file example::

  {
    "files": {
      # An absolute path is left as is
      # a relative path is taken relative to the location of the *script*
      "reference_dir": "/Users/kghose/Data/hg38/",  # Use this if the reference consists of multiple .fa files in a directory
      "reference_file": "/Users/kghose/Data/hg38/hg38.fa.gz",  # Use this if reference is a single multi-fasta file
      "dbfile": "Out/test.h5"    # Genomes database file. Leave out if taking reads from reference
      "output_prefix": "Out/reads", # Output file name prefix
                                    # the reads will be called reads.fq and reads_c.fq if we call for corrupted reads too
      "interleaved": true     # Set to false if you need separate files for each mate of a pair
                              # files will then be named reads_1.fq, reads_c_1.fq and reads_2.fq, reads_c_2.fq
                              # If the reads are actually not paired and you set interleaved to false you will get two
                              # files _1 and _2 and all the data will be in _1 only
      "gzipped": true    # file(s) should be gzipped
    },
    "sample_name": "g0_s0",   # Name of sample
    "rng": {
      "master_seed": 1
    },
    "chromosomes": [1, [2, 0.7, 0.8]],  # List of chromosomes the reads should be taken from
                                        # Note nested notation for taking reads from only part of the chromosome
                                        # [2, 0.7, 0.8]] -> chrom2, take reads only from fraction 0.7 to 0.8 of chromosome
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

__qname_info__ = """
The qname is interpreted as follows:
(Note that the qname is identical for mates from a pair of reads
This is a requirement of the FASTQ convention)

@read_serial|chrom|copy|ro|pos|rlen|cigar|ro|pos|rlen|cigar

read_serial = The read serial is an increasing number
chrom       = chromosome number the read is from
copy        = copy of the chromosome the read is from.
              (Makes a difference for heterozygous variant bearing reads)
                                                                      --@
ro          = 0 - forward strand read (matched to the ref as is)        |
              1 - reverse strand read (has to be reverse complemented)  |
pos         = start position of read in reference sequence (0 indexed)  |
rlen        = length of read                                            |
cigar       = correct CIGAR string                                      |
                                               repeated for each mate --@
"""

import json
import os
import time
import io
import gzip

import numpy as np
import click

import mitty.lib
import mitty.lib.reads as lib_reads
import mitty.lib.io as mio
import mitty.lib.variants as vr

import logging
logger = logging.getLogger(__name__)


class ReadSimulator:
  """A convenience class that wraps the parameters and settings for a read simulation"""
  def __init__(self, base_dir, params, ref_file=None, db_file=None, out_prefix=None):
    """Create a read simulator object

    :param base_dir: the directory with respect to which relative file paths will be resolved
    :param params: dict loaded from json file
    """

    fname_prefix = out_prefix or mitty.lib.rpath(base_dir, params['files']['output_prefix'])
    if not os.path.exists(os.path.abspath(os.path.dirname(fname_prefix))):
      os.makedirs(os.path.dirname(fname_prefix))

    self.ref = mio.Fasta(multi_fasta=ref_file or mitty.lib.rpath(base_dir, params['files'].get('reference_file', None)),
                         multi_dir=mitty.lib.rpath(base_dir, params['files'].get('reference_dir', None)),
                         persistent=True)

    self.sample_name = params.get('sample_name', None)
    if 'dbfile' in params['files'] or db_file is not None:
      pop_db_name = db_file or mitty.lib.rpath(base_dir, params['files']['dbfile'])
      self.pop = vr.Population(fname=pop_db_name, mode='r', in_memory=False)
    else:
      self.pop = None
      logger.debug('Taking reads from reference')

    master_seed = int(params['rng']['master_seed'])
    assert 0 < master_seed < mitty.lib.SEED_MAX
    self.seed_rng = np.random.RandomState(seed=master_seed)

    # In the parameter file we only need to specify a region with the chromosome only if we are not taking reads from
    # the entire chromosome. This is for convenience. Here we have to unpack this into a consistent format
    # We make start_f < stop_f regardless of the order in the parameter file
    self.chromosomes = [ch[0] if type(ch) == list else ch for ch in params['chromosomes']]
    self.chromosome_regions = {
      ch: {'chrom': chr_par[0], 'start_f': min(chr_par[1:]), 'stop_f': max(chr_par[1:])} if type(chr_par) == list else
          {'chrom': chr_par, 'start_f': 0.0, 'stop_f': 1.0} for ch, chr_par in zip(self.chromosomes, params['chromosomes'])}

    self.coverage = float(params['coverage'])

    total_blocks_to_do = self.coverage / float(params['coverage_per_block'])

    chrom_meta = self.ref.get_seq_metadata()
    self.sum_of_chromosome_lengths = float(sum([chrom_meta[c - 1]['seq_len'] for c in self.chromosomes]))
    self.blocks_for_chromosome = {c: int(max(1, round(total_blocks_to_do * (self.chromosome_regions[c]['stop_f'] - self.chromosome_regions[c]['start_f']) * chrom_meta[c - 1]['seq_len'] / self.sum_of_chromosome_lengths)))
                                  for c in self.chromosomes}

    self.read_model = mitty.lib.load_reads_plugin(params['read_model']).Model(**params['model_params'])
    self.corrupt_reads = bool(params['corrupt'])

    self.variants_only = params.get('variants_only', None)
    self.variant_window = int(params.get('variant_window', 200)) if self.variants_only else None

    fname_suffix, open_fun = ('.fq.gz', gzip.open) if params['files'].get('gzipped', False) else ('.fq', open)
    if params['files'].get('interleaved', True):
      self.fastq_fp = [open_fun(fname_prefix + fname_suffix, 'w')] * 2
      self.fastq_c_fp = [open_fun(fname_prefix + '_c' + fname_suffix, 'w')] * 2 if self.corrupt_reads else [None, None]
    else:  # Need two files separately
      self.fastq_fp = [open_fun(fname_prefix + '_1' + fname_suffix, 'w'), open_fun(fname_prefix + '_2' + fname_suffix, 'w')]
      self.fastq_c_fp = [open_fun(fname_prefix + '_c_1' + fname_suffix, 'w'), open_fun(fname_prefix + '_c_2' + fname_suffix, 'w')] if self.corrupt_reads else [None, None]

    self.templates_written = 0
    self.reads_generated = 0
    self.bases_covered = 0

  def get_total_blocks_to_do(self):
    return sum(self.blocks_for_chromosome.values()) * 2  # Two copies for each chromosome

  def get_chromosome_list(self):
    return self.chromosomes

  def get_blocks_to_do(self, chrom):
    return self.blocks_for_chromosome

  def generate_and_save_reads(self, chrom, cpy):
    """Grab the appropriate seq, apply variants as needed, generate reads, roll cigars and then save."""
    if self.pop is not None:
      ml = self.pop.get_variant_master_list(chrom=chrom)
      v_index = self.pop.get_sample_variant_index_for_chromosome(chrom=chrom, sample_name=self.sample_name)
    else:
      ml, v_index = vr.VariantList(), []  # Need a dummy variant list for nulls

    seq, variant_waypoints, var_locs_alt_coords = lib_reads.expand_sequence(self.ref[chrom]['seq'], ml, v_index, cpy)
    seq_c = mitty.lib.string.translate(seq, mitty.lib.DNA_complement)
    coverage_per_block = 0.5 * self.coverage / self.blocks_for_chromosome[chrom]
    for blk in range(self.blocks_for_chromosome[chrom]):
      reads, paired = generate_reads(seq, seq_c, var_locs_alt_coords,
                                     self.variant_window, self.read_model,
                                     coverage_per_block,
                                     self.corrupt_reads,
                                     self.seed_rng,
                                     start_f=self.chromosome_regions[chrom]['start_f'],
                                     stop_f=self.chromosome_regions[chrom]['stop_f'],
                                     variants_only=self.variants_only)
      pos, cigars = lib_reads.roll_cigars(variant_waypoints, reads)
      template_count, bases_covered = write_reads_to_file(
        self.fastq_fp[0], self.fastq_fp[1],
        self.fastq_c_fp[0], self.fastq_c_fp[1],
        reads, paired, pos, cigars,
        chrom, cpy,
        self.templates_written)
      self.templates_written += template_count
      self.reads_generated += reads.shape[0]
      self.bases_covered += bases_covered
      yield

  def get_templates_written(self):
    return self.templates_written

  def get_read_count(self):
    return self.reads_generated

  def get_coverage_done(self):
    return float(self.bases_covered / self.sum_of_chromosome_lengths)
    # This is approximate, since sample will have different length than reference, reference has 'N's, but good enough


def generate_reads(seq, seq_c, var_locs_alt_coords, variant_window,
                   read_model, coverage, corrupt, seed_rng,
                   start_f=0.0, stop_f=1.0,
                   variants_only=False):
  """Wrapper around read function to handle both regular reads as well as reads restricted to around variants

  :param seq:      forward sequence
  :param seq_c:    complement sequence
  :param var_locs_alt_coords: as returned by expand_sequence
  :param variant_window: how many bases before and after variant should we include
  :param read_model: read model object
  :param coverage: real number indicating coverage needed for this run of the simulator
  :param corrupt:  T/F generate corrupted reads too or not
  :param seed_rng:    rng for seed generation
  :param start_f: start fraction for chromosome
  :param stop_f:  stop fraction for chromosome
  :param variants_only: set True if we want reads only from variant regions
  :return:
  """
  start_base = int(len(seq) * start_f)
  stop_base = int(len(seq) * stop_f)
  if variants_only:
    rds, paired = read_model.get_zero_reads()
    reads = [rds]
    for v in var_locs_alt_coords:  # [1:-1]:
      start, stop = max(v - variant_window, 0), min(v + variant_window, len(seq))
      if start > stop_base or stop < start_base: continue
      # v is the pos of the variant in sequence coordinates (rather than ref coordinates)
      these_reads, paired = read_model.get_reads(seq, seq_c,
                                                 start_base=start, end_base=stop,
                                                 coverage=coverage,
                                                 corrupt=corrupt,
                                                 seed=seed_rng.randint(0, mitty.lib.SEED_MAX))
      reads += [these_reads]
    reads = np.concatenate(reads)
  else:
    reads, paired = read_model.get_reads(seq, seq_c,
                                         start_base=start_base, end_base=stop_base,
                                         coverage=coverage,
                                         corrupt=corrupt,
                                         seed=seed_rng.randint(0, mitty.lib.SEED_MAX))
  return reads, paired


def write_reads_to_file(fastq_fp_1, fastq_fp_2,
                        fastq_c_fp_1, fastq_c_fp_2,
                        reads, paired, pos, cigars,
                        chrom, cpy,
                        first_serial_no):
  """
  :param fastq_fp_1:   file pointer to perfect reads file
  :param fastq_fp_2:   file pointer to perfect reads file 2. Same as 1 if interleaving. None if not paired
  :param fastq_c_fp_1: file pointer to corrupted reads file. None if no corrupted reads being written
  :param fastq_c_fp_2: file pointer to corrupted reads file 2. Same as 1 if interleaving. None if no corrupted reads being written. None if not paired
  :param reads:        reads recarray
  :param paired:       bool, are reads paired
  :param pos:          list of POS values for the reads
  :param cigars:       list of cigar values for the reads
  :param chrom:           chromosome number
  :param cpy:           chromosome copy
  :param first_serial_no: serial number of first template in this batch
  :return: next_serial_no: the serial number the next batch should start at

  qname format:

  'read_serial|chrom|copy|ro|pos|rlen|cigar|ro|pos|rlen|cigar'

  ro = readorder => 0 if forward seq, 1 if rev complement

  """
  bases_covered = 0
  ro, pr_seq, cr_seq, phred = reads['read_order'], reads['perfect_reads'], reads['corrupt_reads'], reads['phred']

  if paired:
    for n in xrange(0, reads.shape[0], 2):
      l1, l2 = len(pr_seq[n]), len(pr_seq[n + 1])
      qname = '{:d}|{:d}|{:d}|{:d}|{:d}|{:d}|{:s}|{:d}|{:d}|{:d}|{:s}'.\
        format(first_serial_no + n/2, chrom, cpy, ro[n], pos[n], l1, cigars[n], ro[n + 1], pos[n + 1], l2, cigars[n + 1])
      fastq_fp_1.write('@' + qname + '/1\n' + pr_seq[n] + '\n+\n' + '~' * l1 + '\n')
      fastq_fp_2.write('@' + qname + '/2\n' + pr_seq[n + 1] + '\n+\n' + '~' * l2 + '\n')
      bases_covered += l1 + l2
      if fastq_c_fp_1 is not None:
        fastq_c_fp_1.write('@' + qname + '/1\n' + cr_seq[n] + '\n+\n' + phred[n] + '\n')
        fastq_c_fp_2.write('@' + qname + '/2\n' + cr_seq[n + 1] + '\n+\n' + phred[n + 1] + '\n')
    template_count = reads.shape[0] / 2
  else:
    for n in xrange(0, reads.shape[0]):
      l1 = len(pr_seq[n])
      qname = '{:d}|{:d}|{:d}|{:d}|{:d}|{:d}|{:s}'.format(first_serial_no + n, chrom, cpy, ro[n], pos[n], l1, cigars[n])
      fastq_fp_1.write('@' + qname + '\n' + pr_seq[n] + '\n+\n' + '~' * l1 + '\n')
      bases_covered += l1
      if fastq_c_fp_1 is not None:
        fastq_c_fp_1.write('@' + qname + '\n' + cr_seq[n] + '\n+\n' + phred[n] + '\n')
    template_count = reads.shape[0]
  return template_count, bases_covered


def print_qname(ctx, param, value):
  if not value or ctx.resilient_parsing:
    return
  click.echo(__qname_info__)
  ctx.exit()


@click.group()
@click.version_option()
@click.option('--qname', is_flag=True, callback=print_qname, expose_value=False, is_eager=True, help='Print documentation for information encoded in qname')
def cli():
  """Mitty read simulator"""
  pass


@cli.command()
@click.argument('param_fname', type=click.Path(exists=True))
@click.option('--ref', type=click.Path(exists=True), help="Use this path for reference file. Over-rides entry in parameter file")
@click.option('--db', type=click.Path(exists=True), help="Use this path for genome DB file. Over-rides entry in parameter file")
@click.option('--out-prefix', type=click.Path(), help="Use this path for output file prefix. Over-rides entry in parameter file")
@click.option('-v', count=True, help='Verbosity level')
@click.option('-p', is_flag=True, help='Show progress bar')
def generate(param_fname, ref, db, out_prefix, v, p):
  """Generate reads (fastq) given a parameter file"""
  level = logging.DEBUG if v > 1 else logging.WARNING
  logging.basicConfig(level=level)
  if v == 1:
    logger.setLevel(logging.DEBUG)

  base_dir = os.path.dirname(param_fname)     # Other files will be with respect to this
  params = json.load(open(param_fname, 'r'))

  simulation = ReadSimulator(base_dir, params, ref_file=ref, db_file=db, out_prefix=out_prefix)

  t0 = time.time()
  with click.progressbar(length=simulation.get_total_blocks_to_do(), label='Generating reads', file=None if p else io.BytesIO()) as bar:
    for chrom in simulation.get_chromosome_list():
      for cpy in [0, 1]:
        for _ in simulation.generate_and_save_reads(chrom, cpy):
          bar.update(1)
  t1 = time.time()
  logger.debug('Took {:f}s to write {:d} reads ({:f} coverage)'.format(t1 - t0, simulation.get_read_count(), simulation.get_coverage_done()))


@cli.group()
def show():
  """Show various help pages"""
  pass


@show.command()
def parameters():
  """Program parameter .json"""
  print(__param__)


@show.command('model-list')
def model_list():
  """Print list of available read models"""
  print('\nAvailable models\n----------------')
  for name, mod_name in mitty.lib.discover_all_reads_plugins():
    print('- {:s} ({:s})'.format(name, mod_name))


def explain_all_read_models(ctx, param, value):
  if not value or ctx.resilient_parsing:
    return
  for name, mod_name in mitty.lib.discover_all_reads_plugins():
    explain_read_model(name)
  ctx.exit()


@show.command('read-model')
@click.argument('name')
@click.option('--all', is_flag=True, callback=explain_all_read_models, expose_value=False, is_eager=True, help='Print documentation for all models')
def model(name):
  """Model .json snippets, 'all' for all models"""
  explain_read_model(name)


def explain_read_model(name):
  try:
    mod = mitty.lib.load_reads_plugin(name)
  except ImportError:
    print('No model named {:s}'.format(name))
    return
  try:
    print('\n---- ' + name + ' (' + mod.__name__ + ') ----')
    print(mod._description)
    print(mitty.lib.model_init_signature_string(mod.Model.__init__))
  except AttributeError:
    print('No help for model "{:s}" available'.format(name))
  return


if __name__ == '__main__':
  cli()
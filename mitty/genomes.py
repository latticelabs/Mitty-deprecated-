#!python
from mitty.version import __version__

__cmd__ = """genomes ({:s}): Generate simulated genomes from a simulation parameter file.

Commandline::

  Usage:
    genomes generate <pfile>  [-v|-V] [-p]
    genomes dryrun <pfile>
    genomes inspect <dbfile>
    genomes inspect sfs <chrom> <dbfile>
    genomes write vcf <dbfile> <out_prefix>  [--sample_name=SN] [-v|-V]
    genomes explain (parameters|(variant-model|spectrum-model|population-model) (all|<model_name>))
    genomes list (variant-model|spectrum-model|population-model)

  Options:
    generate                Create a database of genomes by running the models specified
    <pfile>                 Name for parameter file
    dryrun                  Don't run simulation, but print out useful info about it
    inspect                 Give summary information about a genome database
    <dbfile>                Name of genome database file
    sfs                     Print the site frequency spectrum
    <chrom>                 Which chromosome's site frequncy spectrum to print
    write                   Write out genome data from the database file in vcf format
    vcf                     Write out all genomes in one multi-sample vcf file
    <out_prefix>            Written files will have this prefix
    --sample_name=SN        Name of sample. Leave out to write master list
    explain                 Explain the parameters/variant model/population model
    all                     Iterate over and explain all models
    <model_name>            Explain specific model
    list                    List the models
    -v                      Dump log messages
    -V                      Dump detailed log messages
    -p                      Show progress bar""".format(__version__)

__param__ = """Parameter file example::

  {
    # Path notes: an absolute path is left as is. A relative path is taken relative to the parameter file location
    "files": {
      "reference_dir": "/Users/kghose/Data/hg38/",  # If reference is chr1.fa, chr2.fa ... in this directory
      "reference_file": "/Users/kghose/Data/hg38/hg38.fa.gz",  # If reference is a single gz fasta file
      "dbfile": "Out/test.db"  # Output database file
    },
    "rng": {
      "master_seed": 1
    },
    "population_model": {
      "standard": {            # Name of sample chooser (population) model.
        "sample_size": 1,      # population model parameters
      }
    },
    "site_model": {
        "double_exp": {    # Name of model that handles the site frequency spectrum
          "k1": 0.1,       # spectrum model parameters
          "k2": 2.0,
          "p0": 0.001,
          "p1": 0.2,
          "bin_cnt": 30
        }
      }
    "chromosomes": [1, 2]  # Chromosomes to apply the models to
    "variant_models": [    # The list of variant models should come under this key
      {
        "snp": {           # name of the model.
          "p": 0.01        # Parameters required by the model
        }
      },
      {                    # We can chain as many models as we wish. We can also repeat models.
        "delete" : {
          "p": 0.01
        }
      }
    ]
  }"""
__doc__ = __cmd__ + __param__
# We split this up so that when we print the help, we can print it as separate pages. It got cluttered when printed all
# together. We want this to be the docstring because that is a nice format for our generated documentation

import json
import os
import time
import io

import click
import numpy as np

import mitty.lib
import mitty.lib.util as mutil
import mitty.lib.io as mio
import mitty.lib.variants as vr

import logging
logger = logging.getLogger(__name__)


class PopulationSimulator:
  """A convenience class that wraps the parameters and settings for a population simulation"""
  def __init__(self, base_dir, params):
    """Create a genome simulation object

    :param base_dir: the directory with respect to which relative file paths will be resolved
    :param params: dict loaded from json file
    """
    pop_db_name = mitty.lib.rpath(base_dir, params['files']['dbfile'])
    if os.path.exists(pop_db_name):
      logger.warning('Removed old simulation file')
      os.remove(pop_db_name)
    if not os.path.exists(os.path.dirname(pop_db_name)):
      logger.warning('Creating output directory {:s}'.format(pop_db_name))
      os.makedirs(os.path.dirname(pop_db_name))

    self.ref = mio.Fasta(multi_fasta=mitty.lib.rpath(base_dir, params['files'].get('reference_file', None)),
                    multi_dir=mitty.lib.rpath(base_dir, params['files'].get('reference_dir', None)))  # TODO: Ability to switch off persistence flag
    master_seed = int(params['rng']['master_seed'])
    assert 0 < master_seed < mitty.lib.SEED_MAX

    self.seed_rng = np.random.RandomState(seed=master_seed)
    self.pop = vr.Population(fname=pop_db_name, genome_metadata=self.ref.get_seq_metadata())

    self.chromosomes = params['chromosomes']
    self.sfs_model = load_site_frequency_model(params.get('site_model', None))
    self.sfs_p, self.sfs_f = self.sfs_model.get_spectrum() if self.sfs_model is not None else (None, None)
    self.variant_models = load_variant_models(self.ref, params['variant_models'])
    self.population_model = load_population_model(params.get('population_model', None), params)

    self.unique_variant_count, self.total_variant_count = 0, 0

  def get_chromosome_list(self):
    return self.chromosomes

  def get_total_blocks_to_do(self):
    return len(self.chromosomes) * self.population_model.get_sample_count_estimate()

  def generate_and_save_samples(self, chrom):
    ml = vr.VariantList()
    for m in self.variant_models:
      ml.add(*m.get_variants(ref=self.ref[chrom]['seq'], chrom=chrom,
                             p=self.sfs_p, f=self.sfs_f,
                             seed=self.seed_rng.randint(mutil.SEED_MAX)))
    ml.sort()
    if self.sfs_model is not None: ml.balance_probabilities(*self.sfs_model.get_spectrum())
    self.pop.set_master_list(chrom=chrom, master_list=ml)
    self.unique_variant_count += len(ml)
    for sample_name, this_sample, frac_done in self.population_model.samples(chrom_no=chrom, ml=ml, rng_seed=self.seed_rng.randint(mutil.SEED_MAX)):
      self.pop.add_sample_chromosome(chrom=chrom, sample_name=sample_name, indexes=this_sample)
      self.total_variant_count += len(this_sample)
      yield


def load_site_frequency_model(sfs_model_json):
  if sfs_model_json is None:
    return None
  k, v = sfs_model_json.items()[0]
  return mitty.lib.load_sfs_plugin(k).Model(**v)


def load_variant_models(ref, model_param_json):
  """Given a json snippet corresponding to models and their parameters, load them"""
  return [mitty.lib.load_variant_plugin(k).Model(ref=ref, **v)
          for model_json in model_param_json
          for k, v in model_json.iteritems()]  # There really is only one key (the model name) and the value is the
                                               # parameter list


def load_population_model(pop_model_json, params={}):
  """Given a json snippet corresponding to the population model load it. If None, return the standard model

  :param pop_model_json: json snippet corresponding to the population model
  :param params:  the entire parameter json. For backward compatibility. If we have no pop model defined, as was normal
                  for Mitty v < 1.2.0, we default to the standard model (built-in). The single parameter, sample_size,
                  was defined in the main body of the parameter json, so we need that here.
  """
  if pop_model_json is None:
    pop_model_json = {'standard': {'sample_size': params['sample_size']}}  # params is assumed to have the key 'sample_size'
  k, v = pop_model_json.items()[0]
  return mitty.lib.load_pop_model_plugin(k).Model(**v)


@click.group()
@click.version_option()
def cli():
  """Mitty genomes simulator"""
  logging.basicConfig(level=logging.INFO)


@cli.command()
@click.argument('param_fname', type=click.Path(exists=True))
@click.option('--dry-run', is_flag=True, help='Print useful information about simulation, but dont run')
@click.option('-v', count=True, help='Verbosity level')
@click.option('-p', is_flag=True, help='Show progress bar')
def generate(param_fname, dry_run, v, p):
  """Generate population of genomes"""
  level = logging.DEBUG if v > 1 else logging.WARNING
  logging.basicConfig(level=level)
  if v == 1:
    logger.setLevel(logging.DEBUG)

  base_dir = os.path.dirname(param_fname)     # Other files will be with respect to this
  params = json.load(open(param_fname, 'r'))

  if dry_run:
    do_dry_run(params)
    return

  simulation = PopulationSimulator(base_dir, params)
  t0 = time.time()
  with click.progressbar(length=simulation.get_total_blocks_to_do(), label='Generating genomes', file=None if p else io.BytesIO()) as bar:
    for chrom in simulation.get_chromosome_list():
      for _ in simulation.generate_and_save_samples(chrom):
        bar.update(1)
  t1 = time.time()
  logger.debug('Took {:f}s'.format(t1 - t0))
  logger.debug('{:d} unique variants, {:d} variants in samples'.format(simulation.unique_variant_count, simulation.total_variant_count))


def do_dry_run(params):
  """Print useful info about simulation"""
  sfs_model = load_site_frequency_model(params.get('site_model', None))
  if sfs_model:
    print('Site frequency model')
    print(sfs_model)
    p, f = sfs_model.get_spectrum()
    if len(p):  # Don't print this for a dud spectrum
      print('1/sum(p_i * f_i) = {:2.1f}'.format(1./(p*f).sum()))
    x, y = mutil.growth_curve_from_sfs(p, f)
    print('\nGrowth curve:\n')
    print('N\tv')
    for _x, _y in zip(x, y):
      print('{:d}\t{:f}'.format(_x, _y))
  else:
    print('No site model')


@cli.group('genome-file')
def g_file():
  """Inspect or operate on existing genome file."""
  pass


@g_file.command('write-vcf')
@click.argument('dbfile', type=click.Path(exists=True))
@click.argument('vcfgz', type=click.Path())
@click.option('--sample-name', help='Name of sample. Omit to write master list')
def write_vcf(dbfile, vcfgz, sample_name):
  """Write sample/master list to VCF"""
  pop = vr.Population(fname=dbfile)
  if sample_name in pop.get_sample_names() or sample_name is None:
    mio.write_single_sample_to_vcf(pop=pop, sample_name=sample_name, out_fname=vcfgz)
  else:
    logger.warning('Sample name {:s} not in population'.format(sample_name))
  # if sample_name is none, mio.write_single_sample_to_vcf will write master list


@g_file.command('summary')
@click.argument('dbfile', type=click.Path(exists=True))
def summary(dbfile):
  """Print some useful information about the database

  :param cmd_args: parsed arguments
  """
  print(vr.Population(dbfile))


@g_file.command('sfs')
@click.argument('dbfile', type=click.Path(exists=True))
@click.argument('chrom', type=int)
def print_sfs(dbfile, chrom):
  """Print site frequency spectrum for chrom in file"""
  pop = vr.Population(fname=dbfile)
  ml = pop.get_master_list(chrom)
  print('Site frequency spectrum for chrom {:d}'.format(chrom))
  print(ml)


@cli.group()
def show():
  """Various help pages"""
  pass

@show.command()
def parameters():
  """Program parameter .json"""
  print(__param__)


@show.command('model-list')
def model_list():
  """Print list of models"""
  discoverer = [
    ('variant', mitty.lib.discover_all_variant_plugins),
    ('spectrum', mitty.lib.discover_all_sfs_plugins),
    ('population', mitty.lib.discover_all_pop_plugins)
  ]
  for mod_type, disco in discoverer:
    print('--------------------------------\nAvailable {} models\n--------------------------------'.format(mod_type))
    for name, mod_name in disco():
      print('- {:s} ({:s})'.format(name, mod_name))

@show.command('variant-model')
@click.argument('name')
def explain_variant_models(name):
  """Explain variant models. 'all' for all models"""
  if name == 'all':
    explain_all_models('variant-model')
  else:
    explain_model(name, 'variant-model')


@show.command('spectrum-model')
@click.argument('name')
def explain_spectrum_models(name):
  """Explain spectrum models. 'all' for all models"""
  if name == 'all':
    explain_all_models('spectrum-model')
  else:
    explain_model(name, 'spectrum-model')


@show.command('population-model')
@click.argument('name')
def explain_population_models(name):
  """Explain population models. 'all' for all models"""
  if name == 'all':
    explain_all_models('population-model')
  else:
    explain_model(name, 'population-model')


def explain_all_models(kind):
  discoverer = {
    'variant-model': mitty.lib.discover_all_variant_plugins,
    'spectrum-model': mitty.lib.discover_all_sfs_plugins,
    'population-model': mitty.lib.discover_all_pop_plugins
  }
  for name, mod_name in discoverer[kind]():
    explain_model(name, kind)


def explain_model(name, kind):
  """Load the given model of given kind and print description, parameter example and defaults

  :param name: name of model
  :param kind: 'variant', 'sfs' or 'pop'
  """
  loader = {
    'variant-model': mitty.lib.load_variant_plugin,
    'spectrum-model': mitty.lib.load_sfs_plugin,
    'population-model': mitty.lib.load_pop_model_plugin
  }
  try:
    mod = loader[kind](name)
  except ImportError as e:
    print('{0}: {1}'.format(name, e))
    print('Problem with loading {} model'.format(kind))
    return
  try:
    print('\n---- ' + name + ' (' + mod.__name__ + ') ----')
    print(mod._description)
    print(mitty.lib.model_init_signature_string(mod.Model.__init__))
  except AttributeError:
    print('No help for model "{:s}" available'.format(name))


if __name__ == "__main__":
  cli()
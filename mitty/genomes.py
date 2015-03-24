#!python
__cmd__ = """Generate simulated genomes from a simulation parameter file.

Commandline::

  Usage:
    genomes generate --pfile=PFILE  [-v|-V]
    genomes write (vcf|vcfs) --dbfile=DBFILE  (<serial>)... [-v|-V]
    genomes explain (parameters|(variantmodel|populationmodel) <model_name>)
    genomes list (variantmodels|populationmodels)

  Options:
    generate                Create a database of genomes by running the models specified
    --pfile=PFILE           Name for parameter file
    -v                      Dump log messages
    -V                      Dump detailed log messages
    write                   Write out genome data from the database file in vcf format
    vcf                     Write out all genomes in one multi-sample vcf file
    vcfs                    Write out the genomes in separate single sample vcf files
    --dbfile=DBFILE         Name of genome database file
    <serial>                Serial number of sample
    explain                 Explain the parameters/variant model/population model
    list                    List the models
"""

__param__ = """Parameter file example::

  {
    # Path notes: an absolute path is left as is. A relative path is taken relative to the parameter file location
    "files": {
      "reference_dir": "/Users/kghose/Data/hg38/",  # Use this if the reference consists of multiple .fa files in a directory
      "reference_file": "/Users/kghose/Data/hg38/hg38.fa.gz",  # Use this if reference is a single multi-fasta file
      "dbfile": "Out/test.db"  # Output database file
    },
    "rng": {
      "master_seed": 1
    },
    "sample_size": 1,            # How many samples to generate
    "site_model": {
        "double_exp": {   # Name of model that handles the site frequency spectrum
          "k1": 0.1,     # Population model parameters
          "k2": 2.0,
          "p0": 0.001,
          "p1": 0.2,
          "bin_cnt": 30
        }
      }
    "chromosomes": [1, 2]        # Chromosomes to apply the models to
    "variant_models": [          # The list of variant models should come under this key
      {
        "snp": {                 # name of the model. To get a list of plugin names type "denovo models"
          "p": 0.01              # Parameters required by the model
        }
      },
      {                          # We can chain as many models as we wish
        "delete" : {             # We can repeat models if we want
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
from collections import deque
import time

import docopt
import numpy as np

import mitty.lib
import mitty.lib.util as mutil
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

  if cmd_args['generate']:
    generate(cmd_args)
  elif cmd_args['write']:
    write(cmd_args)
  elif cmd_args['explain']:
    explain(cmd_args)
  elif cmd_args['list']:
    print_list(cmd_args)


def generate(cmd_args):
  """Generate genomes based on the simulation parameter file"""

  base_dir = os.path.dirname(cmd_args['--pfile'])     # Other files will be with respect to this
  params = json.load(open(cmd_args['--pfile'], 'r'))

  #conn = db.connect(mitty.lib.rpath(base_dir, params['files']['dbfile']))  # We save our data to this db

  ref_genome = mio.Fasta(multi_fasta=params['files'].get('reference_file', None),
                         multi_dir=params['files'].get('reference_dir', None))  # TODO: Ability to switch off persistence flag
  master_seed = int(params['rng']['master_seed'])
  assert 0 < master_seed < mitty.lib.SEED_MAX

  pop_size = int(params['population_size'])
  sample_size = int(params['sample_size'])





  generations = deque()
  ml = {}

  seeds = mitty.lib.get_seeds(master_seed, n_generations)

  g1 = gen.create_initial_generation(init_pop_size)
  models = load_variant_models(ref_genome, params.get('initial_pop_variant_models', []))
  add_denovo_variants_from_models_to_population(g1, models, ml, seeds[0])
  generations.append(g1)

  models = load_variant_models(ref_genome, params.get('steady_state_variant_models', []))
  for n, seed in enumerate(seeds[1:]):
    logger.debug('Running generation {:d} ({:d})'.format(n, len(g1)))
    logger.debug(g1)
    g1 = gen.create_next_generation(g1, fitness, mater, breeder, culler)
    gen.genomify(g1, crosser)
    add_denovo_variants_from_models_to_population(g1, models, ml, seed)
    # generations.append(g1)
    # if n > n_ancestors:
    #   _ = generations.popleft()  # Get rid of earliest ancestors to save space

  #db.save_variant_master_list(conn, ml)


def run_simulations(pop_db_name, ref, sfs_model, variant_models=[], chromosomes=[], sample_size=1, master_seed=2):
  seed_rng = np.random.RandomState(seed=master_seed)
  conn = mdb.connect(db_name=pop_db_name)
  for ch in chromosomes:
    ml = vr.VariantList()
    for m in variant_models:
      ml.add(*m.get_variants(ref[ch], ch, *sfs_model.get_spectrum(), seed=seed_rng.randint(mutil.SEED_MAX)))
    ml.sort()
    ml.balance_probabilities(*sfs_model.get_spectrum())
    mdb.save_master_list(conn, ch, ml)
    rng = np.random.RandomState(seed_rng.randint(mutil.SEED_MAX))
    for n in range(sample_size):
      mdb.save_sample(conn, 0, n, ch, ml.generate_chromosome(rng))
  conn.close()


def load_variant_models(ref, model_param_json):
  """Given a list of models and parameters load the relevant modules and store the parameters as tuples"""
  return [mitty.lib.load_variant_plugin(k).Model(ref=ref, **v)
          for model_json in model_param_json
          for k, v in model_json.iteritems()]  # There really is only one key (the model name) and the value is the
                                               # parameter list


def add_denovo_variants_from_models_to_genome(g1, models, ml, master_seed=1):
  """Given an original genome add any variants that come off the variant_generator. g1 is modified in place

  :param dict g1: genome
  :param variant_generators: list of variant generator iterators
  :param ml: master list of variants
  :returns: g1 is modified in place
  """
  for model, seed in zip(models, mitty.lib.get_seeds(master_seed, len(models))):
    dnv = model.variants(seed=seed)
    # dnv is a dict with keys as chrom numbers and items as list of lists (pos, ref, alt, gt) from the variant plugin
    for chrom, pvd in dnv.iteritems():  # pvd is a list of lists (pos, ref, alt, gt) from the variant plugin
      dnv = vr.create_gtv_iterable(pvd[0], pvd[1], pvd[2], pvd[3])
      if chrom not in g1:  # Need to make a new chromosome for this.
        g1[chrom] = vr.Chromosome()
      if chrom not in ml:  # Need to make a new chromosome for this.
        ml[chrom] = {}
      vr.add_denovo_variants_to_chromosome(g1[chrom], dnv, ml[chrom])


def add_denovo_variants_from_models_to_population(p, models, ml, master_seed=1):
  """Given a population (list of genomes) and a model list, repeatedly apply the models to each genome

  :param p: list of gen.Sample from the same generation
  :param mdl_list: list of models
  :param ref: reference genome
  :param ml: master list of variants
  :param master_seed: this will drive seeds for this part of the simulation
  :returns: modifies all genomes in p, in place
  """
  for g, seed in zip(p, mitty.lib.get_seeds(master_seed, len(p))):
    add_denovo_variants_from_models_to_genome(g.genome, models, ml, seed)


def write(cmd_args):
  pass


def explain(cmd_args):
  if cmd_args['parameters']:
    print(__param__)
  elif cmd_args['variantmodel']:
    explain_variant_model(cmd_args['<model_name>'])
  elif cmd_args['populationmodel']:
    explain_population_model(cmd_args['<model_name>'])


def explain_variant_model(name):
  try:
    mod = mitty.lib.load_variant_plugin(name)
  except ImportError as e:
    print('{0}: {1}'.format(name, e))
    print('Problem with loading model')
    return
  try:
    print('\n---- ' + name + ' (' + mod.__name__ + ') ----')
    print(mod._description)
  except AttributeError:
    print('No help for model "{:s}" available'.format(name))


def explain_population_model(name):
  pass


def print_list(cmd_args):
  if cmd_args['variantmodels']:
    print_variant_model_list()
  elif cmd_args['populationmodels']:
    print_population_model_list()


def print_variant_model_list():
  print('\nAvailable variant models\n----------------')
  for name, mod_name in mitty.lib.discover_all_variant_plugins():
    print('- {:s} ({:s})\n'.format(name, mod_name))


def print_population_model_list():
  pass

  #   print_model_list()
  #   exit(0)
  # if cmd_args['explain']:
  #   if '<model_name>' in cmd_args:
  #     explain_model(cmd_args['<model_name>'])
  #   exit(0)
  #
  # level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  # logging.basicConfig(level=level)
  #
  # import os
  # base_dir = os.path.dirname(cmd_args['--pfile'])     # Other files will be with respect to this
  # params = json.load(open(cmd_args['--pfile'], 'r'))
  # ref_genome = FastaGenome(seq_dir=mitty.lib.rpath(base_dir, params['files']['genome']), persist=True)
  # vcf_file_name = mitty.lib.rpath(base_dir, params['files']['output vcf'])
  # models = load_variant_model_list(params['denovo_variant_models'])
  # master_seed = params['rng']['master_seed']
  #
  # t0 = time.time()
  # g1 = main(ref=ref_genome, models=models, master_seed=master_seed)
  # t1 = time.time()
  # logger.debug('Computed variants in {:f} s'.format(t1 - t0))
  #
  # t0 = time.time()
  # mitty.lib.io.vcf_save_gz(g1, vcf_file_name)
  # t1 = time.time()
  # logger.debug('Saved VCF file in {:f} s'.format(t1 - t0))
  #
  #

if __name__ == "__main__":
  cli()
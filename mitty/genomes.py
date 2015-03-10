#!python
__cmd__ = """This command line program generates a database of genomes given a simulation parameter file.

Commandline::

  Usage:
    genomes generate --pfile=PFILE  [--continue] [-v]
    genomes write (vcf|vcfs) --dbfile=DBFILE  (<gen> <serial>)... [-v]
    genomes explain (parameters|variantmodel <model_name>|populationmodel <model_name>)
    genomes list (variantmodels|populationmodels)

  Options:
    generate                Create a database of genomes by running the models specified
    --pfile=PFILE           Name for parameter file
    --continue              Continue on from an interrupted population simulation
    -v                      Dump detailed logger messages
    write                   Write out the genomes
    vcf                     Write out genomes in one multi-sample vcf file
    vcfs                    Write out the genomes in separate single sample vcf files
    --dbfile=DBFILE         Name of genome database file
    gen                     Generation of sample
    serial                  Serial number of sample
    explain                 Explain something
    parameters              Pring an example parameter file
    variantmodel            Access the variant model plugins
    <model_name>            The model to explain.
    populationmodel         Access the population model plugins
    list                    List the models
"""

__param__ = """Parameter file example::

  {
    "files": {
      "reference": "/Users/kghose/Data/hg38/",  # An absolute path is left as is
      "dbfile": "Out/test.vcf.gz"          # a relative path is taken relative to the location of the *script*
    },
    "rng": {
      "master_seed": 1
    },
    "generations": 10            # Number of generations to run the model
    "initial_pop_size": 10       # Number of genomes to start with, created denovo
    "generations_to_keep": 1     # How many generations (in addition to the current one) should we store in the database
    "initial_pop_variant_models": [   # The list of variant models should come under this key
      {
        "snp": {                 # name of the model. To get a list of plugin names type "denovo models"
          "chromosome": [1, 2],  # Chromosomes to apply this model to
          "phet": 0.5,           # Parameters required by the model
          "p": 0.01,
          "poisson_rng_seed": 1,
          "base_sub_rng_seed": 2
        }
      },
      {                          # We can chain as many models as we wish
        "delete" : {             # We can repeat models if we want
          "chromosome": [1],
          "phet": 0.5,
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

import mitty.lib
import mitty.lib.io as mio
#from mitty.lib.genome import FastaGenome
import mitty.lib.variation as vr

import logging
logger = logging.getLogger(__name__)


def cli():  # pragma: no cover
  """Serves as entry point for scripts"""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__cmd__, ['-h'])
  else:
    cmd_args = docopt.docopt(__cmd__)  # Version string?

  level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  if cmd_args['generate']:
    generate(cmd_args)
  elif cmd_args['write']:
    write(cmd_args)
  elif cmd_args['explain']:
    explain(cmd_args)
  elif cmd_args['list']:
    print_list(cmd_args)


def generate(cmd_args):
  """Run a simulation based on the parameter file."""

  base_dir = os.path.dirname(cmd_args['--pfile'])     # Other files will be with respect to this
  params = json.load(open(cmd_args['--pfile'], 'r'))

  ref_genome = mio.load_multi_fasta_gz(mitty.lib.rpath(base_dir, params['files']['reference']))
  master_seed = int(params['rng']['master_seed'])
  assert 0 < master_seed < mitty.lib.SEED_MAX

  n_generations = int(params['generations'])
  n_ancestors = int(params['generations_to_keep'])

  generations = deque()
  ml = {}

  seeds = mitty.lib.get_seeds(master_seed, n_generations)

  g1 = [{} for _ in range(int(params['initial_pop_size']))]
  mdl_list = load_variant_model_list(params.get('initial_pop_variant_models', []))
  add_denovo_variants_from_models_to_population(g1, mdl_list, ref_genome, ml, seeds[0])
  generations.append(g1)

  mdl_list = load_variant_model_list(params.get('steady_state_variant_models', []))
  for n, seed in enumerate(seeds[1:]):
    g1 = create_next_generation(g1)
    add_denovo_variants_from_models_to_population(g1, mdl_list, ref_genome, ml, seed)
    generations.append(g1)
    if n > n_ancestors:
      _ = generations.popleft()  # Get rid of earliest ancestors to save space


def load_variant_model_list(model_param_json):
  """Given a list of models and parameters load the relevant modules and store the parameters as tuples"""
  return [{"model": mitty.lib.load_variant_plugin(k), "params": v}
          for model_json in model_param_json
          for k, v in model_json.iteritems()]  # There really is only one key (the model name) and the value is the
                                               # parameter list


def get_variant_generator_list_from_model_list(mdl_list, ref, master_seed=1):
  """Initialize the variant generators from the list of models (and parameters) while passing a master seed."""
  assert 0 < master_seed < mitty.lib.SEED_MAX
  return [mdl['model'].variant_generator(ref=ref, master_seed=seed, **mdl['params'])
          for mdl, seed in zip(mdl_list, mitty.lib.get_seeds(master_seed, len(mdl_list)))]


def add_denovo_variants_from_models_to_genome(g1, variant_generators, ml):
  """Given an original genome add any variants that come off the variant_generator. g1 is modified in place

  :param dict g1: genome
  :param variant_generators: list of variant generator iterators
  :param ml: master list of variants
  :returns: g1 is modified in place
  """
  for vg in variant_generators:
    for delta_g in vg:  # delta_g is a dictionary of proposed variants, indexed by chromosome number (1,2....)
      for chrom, pvd in delta_g.iteritems():  # pvd is a list of lists (pos, ref, alt, gt) from the variant plugin
        dnv = vr.create_sample_iterable(pvd[0], pvd[1], pvd[2], pvd[3])
        if chrom not in g1:  # Need to make a new chromosome for this.
          g1[chrom] = vr.Sample()
        s = g1.get(chrom)
        vr.add_denovo_variants_to_sample(s, dnv, ml)


def add_denovo_variants_from_models_to_population(p, mdl_list, ref, ml, master_seed=1):
  """Given a population (list of genomes) and a model list, repeatedly apply the models to each genome

  :param p: list of genomes from the same generation
  :param mdl_list: list of models
  :param ref: reference genome
  :param ml: master list of variants
  :param master_seed: this will drive seeds for this part of the simulation
  :returns: modifies all genomes in p, in place
  """
  for g, seed in zip(p, mitty.lib.get_seeds(master_seed, len(p))):
    v_gen = get_variant_generator_list_from_model_list(mdl_list, ref, seed)
    add_denovo_variants_from_models_to_genome(g, v_gen, ml)


def create_next_generation(g0):
  pass


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
  except ImportError:
    print('No model named {:s}'.format(name))
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








def main(ref, models=[], master_seed=1):
  """This does what the old mutate.py script used to do."""
  assert 0 < master_seed < mitty.lib.SEED_MAX
  logger.debug('Reference file {:s}'.format(ref.dir))
  vgl = get_variant_generator_list_from_model_list(models, ref, master_seed=master_seed)
  return merge_variants_from_models(g1={}, variant_generators=vgl)




if __name__ == "__main__":
  cli()
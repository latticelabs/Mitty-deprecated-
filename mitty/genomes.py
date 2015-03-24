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

  pop_db_name = mitty.lib.rpath(base_dir, params['files']['dbfile'])

  ref = mio.Fasta(multi_fasta=mitty.lib.rpath(base_dir, params['files'].get('reference_file', None)),
                  multi_dir=mitty.lib.rpath(base_dir, params['files'].get('reference_dir', None)))  # TODO: Ability to switch off persistence flag
  master_seed = int(params['rng']['master_seed'])
  assert 0 < master_seed < mitty.lib.SEED_MAX

  sample_size = int(params['sample_size'])
  chromosomes = params['chromosomes']

  run_simulations(pop_db_name, ref, sfs_model=load_site_frequency_model(params['site_model']),
                  variant_models=load_variant_models(ref, params['variant_models']),
                  chromosomes=chromosomes, sample_size=sample_size, master_seed=master_seed)


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


def load_site_frequency_model(sfs_model_json):
  k, v = sfs_model_json.items()[0]
  return mitty.lib.load_sfs_plugin(k).Model(**v)


def load_variant_models(ref, model_param_json):
  """Given a json snippet corresponding to models and their parameters, load them"""
  return [mitty.lib.load_variant_plugin(k).Model(ref=ref, **v)
          for model_json in model_param_json
          for k, v in model_json.iteritems()]  # There really is only one key (the model name) and the value is the
                                               # parameter list


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


if __name__ == "__main__":
  cli()
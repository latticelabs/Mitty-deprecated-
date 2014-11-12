#!python
"""This script produces a VCF file simulating a "sample" based on a reference

Commandline::

  Usage:
    denovo --pfile=PFILE  [-v]
    denovo plugins
    denovo explain <plugin>

  Options:
    --pfile=PFILE           Name for parameter file
    -v                      Dump detailed logger messages
    plugins                 List the available denovo plugins
    explain                 Explain details about the indicated plugin
    <plugin>                The plugin to explain. If none, explains the parameter file format


Parameter file example::

  {
    "files": {
      "genome": "/Users/kghose/Data/hg38",  # An absolute path is left as is
      "output vcf": "Out/test.vcf.gz"          # a relative path is taken relative to the location of the *script*
    },
    "rng": {
      "master_seed": 1
    },
    "denovo_variant_models": [    # The list of variant models should come under this key
      {
        "snp": {                 # name of the model. To get a list of plugin names type "denovo plugins"
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
  }

Notes:

1. denovo sorts, compresses and indexes the output vcf, which should be given as a vcf.gz file
"""
__version__ = '1.0.0'

import numpy
import json
import docopt
from lib.genome import FastaGenome
from mitty.lib.variation import *
from mitty.lib import *
import logging
from plugins import putil
logger = logging.getLogger(__name__)


def merge_variants_with_genome(g1, variant_generator):
  """Given an original genome add any variants that come off the variant_generator
  Args:
    g1 (dict): Dict of chromosomes
  """
  g2 = copy_genome(g1)
  for delta_g in variant_generator:
    for chrom, dnv in delta_g.iteritems():
      try:
        g2[chrom] = merge_variants_with_chromosome(g2[chrom], dnv)
      except KeyError:  # This is the first time we are seeing variants on this chromosome
        g2[chrom] = dnv
  return g2


def apply_variant_models_to_genome(g1={}, ref=None, models=[]):
  """Return a copy of g1 with the denovo variants from the models merged in"""
  g2 = None
  for model in models:
    g2 = merge_variants_with_genome(g2 or g1, model["model"].variant_generator(ref, **model["params"]))
  return g2


def create_denovo_genome(ref, models=[]):
  """Simply a convenience function that creates a genome from scratch."""
  return apply_variant_models_to_genome(ref=ref, models=models)


def load_variant_models(model_param_json):
  """Given a list of models and parameters load the models in"""
  return [{"model": putil.load_variant_plugin(k), "params": v}
          for model_json in model_param_json
          for k, v in model_json.iteritems()]  # There really is only one key (the model name) and the value is the
                                               # parameter list


def apply_master_seed(models, master_seed=1):
  assert 0 < master_seed < SEED_MAX
  for model, seed in zip(models, numpy.random.RandomState(seed=master_seed).randint(SEED_MAX, size=len(models))):
    model["params"]["master_seed"] = seed


def print_plugin_list():
  print 'Available plugins'
  for plugin in putil.list_all_variant_plugins():
    print plugin


def explain_plugin(plugin):
  if plugin not in putil.list_all_variant_plugins():
    print 'No such plugin'
  else:
    mod = putil.load_variant_plugin(plugin)
    try:
      print(mod._description)
    except AttributeError:
      print('No help for module available')
  return


def main(ref, vcf_file_name=None, parameters={}, master_seed=1):
  """This does what the old mutate.py script used to do."""
  assert 0 < master_seed < SEED_MAX
  logger.debug('Reference file {:s}'.format(ref.dir))
  models = load_variant_models(parameters['denovo_variant_models'])
  apply_master_seed(models, master_seed)
  g1 = create_denovo_genome(ref=ref, models=models)
  if vcf_file_name is not None:
    vcf_save_gz(g1, vcf_file_name)
  return g1


if __name__ == "__main__": # pragma: no cover
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__, version=__version__)
  if cmd_args['plugins']:
    print_plugin_list()
    exit(0)
  if cmd_args['explain']:
    explain_plugin(cmd_args['<plugin>'])
    exit(0)

  level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  base_dir = os.path.dirname(cmd_args['--pfile'])     # Other files will be with respect to this
  params = json.load(open(cmd_args['--pfile'], 'r'))
  ref_genome = FastaGenome(seq_dir=rpath(base_dir, params['files']['genome']), persist=True)
  vcf_file = params['files']['output vcf']
  vcf_file = rpath(base_dir, vcf_file) if vcf_file else vcf_file
  master_seed = params['rng']['master_seed']
  main(ref=ref_genome, vcf_file_name=vcf_file, parameters=params, master_seed=master_seed)
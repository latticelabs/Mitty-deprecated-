#!python
"""This module implements functions to add denovo variants to a genome. When run as a script it will produce a VCF
file simulating a "sample".

Commandline::

  Usage:
    denovo --pfile=PFILE  [-v]
    denovo models
    denovo explain <model_name>

  Options:
    --pfile=PFILE           Name for parameter file
    -v                      Dump detailed logger messages
    models                  List the available denovo models
    explain                 Explain details about the indicated plugin
    <model_name>            The model to explain. If none, explains the parameter file format


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

import json
import docopt
import mitty.lib
import mitty.lib.io
from mitty.lib.genome import FastaGenome
from mitty.lib.variation import *
import logging
logger = logging.getLogger(__name__)


def copy_missing_chromosomes(g1, g2):
  """Copy any chromosomes found in g2 but not in g1 onto g1. g1 is changed in place

  :param dict g1: genome
  :param dict g2: genome
  :returns: changes g1 in place"""
  missing_chrom = set(g2.keys()) - set(g1.keys())
  for ch in missing_chrom:
    g1[ch] = copy_variant_sequence(g2[ch])


def merge_variants_from_models(g1={}, variant_generators=[]):
  """Given an original genome add any variants that come off the variant_generator

  :param dict g1: genome
  :param list variant_generators: list of variant generators
  :returns: list of Variants
  """
  g2 = {}
  for vg in variant_generators:
    for delta_g in vg:
      for chrom, dnv in delta_g.iteritems():
        g2[chrom] = merge_variants(g2.get(chrom) or g1.get(chrom, []), dnv)
  copy_missing_chromosomes(g2, g1)
  return g2


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


def create_variant_list_from_models(mdl_list, ref, master_seed=1):
  """Given a model list, seed and ref return us a genome. Convenience wrapper function."""
  return merge_variants_from_models(variant_generators=
                                    get_variant_generator_list_from_model_list(mdl_list, ref, master_seed))


def print_model_list():
  print('\nAvailable models\n----------------')
  for name, mod_name in mitty.lib.discover_all_variant_plugins():
    print('- {:s} ({:s})\n'.format(name, mod_name))


def explain_model(name):
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
  return


def main(ref, models=[], master_seed=1):
  """This does what the old mutate.py script used to do."""
  assert 0 < master_seed < mitty.lib.SEED_MAX
  logger.debug('Reference file {:s}'.format(ref.dir))
  vgl = get_variant_generator_list_from_model_list(models, ref, master_seed=master_seed)
  return merge_variants_from_models(g1={}, variant_generators=vgl)


def cli():  # pragma: no cover
  """Serves as entry point for scripts"""
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__, version=__version__)
  if cmd_args['models']:
    print_model_list()
    exit(0)
  if cmd_args['explain']:
    explain_model(cmd_args['<model_name>'])
    exit(0)

  level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  import os
  base_dir = os.path.dirname(cmd_args['--pfile'])     # Other files will be with respect to this
  params = json.load(open(cmd_args['--pfile'], 'r'))
  ref_genome = FastaGenome(seq_dir=mitty.lib.rpath(base_dir, params['files']['genome']), persist=True)
  vcf_file_name = mitty.lib.rpath(base_dir, params['files']['output vcf'])
  models = load_variant_model_list(params['denovo_variant_models'])
  master_seed = params['rng']['master_seed']

  g1 = main(ref=ref_genome, models=models, master_seed=master_seed)
  mitty.lib.io.vcf_save_gz(g1, vcf_file_name)


if __name__ == "__main__":
  cli()
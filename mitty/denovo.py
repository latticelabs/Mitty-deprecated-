"""This module contains the machinery to use variation plugins to generate de novo variants to a given reference

This will gradually replace mutate, since it is effectively a superset of mutate. Mutate has a cooler name tho - maybe
mutate will be the command line script that wraps population functions ...



This module contains functions that generate variants with reference to a reference genome. The module can be called
as a script as well. This is useful for creating test data for MGR algorithms/data formats. The script will output
VCF file(s).

Commandline::

  Usage:
    denovo --wg=WG  --vcf=VCF  --param_file=PFILE  [--master_seed=SEED] [-v]

  Options:
    --wg=WG                 Whole genome reference file
    --vcf=VCF               Output VCF file. mutate expects the file to end in .vcf.gz
                            It saves this file as a sorted and compressed vcf file with the uncompressed version
                            available as .vcf. If the file does not end in .vcf.gz the uncompressed version will
                            be of the form <basename>_srt.vcf
    --paramfile=PFILE       Name for parameter file
    --master_seed=SEED      If this is specified, this generates and passes master seeds to all the plugins.
                            This overrides any individual seeds specified by the parameter file.
    -v                      Dump detailed logger messages

Parameter file example::

  {
    "variant_models": [
      {
        "chromosome": [1],
        "model": "snp",
        "phet": 0.5,
        "p": 0.01,
        "poisson_rng_seed": 1,
        "base_sub_rng_seed": 2
      },
      {
        "chromosome": [1],
        "model": "delete",
        "phet": 0.5,
        "p_del": 0.01,Re
        "lam_del": 10,
        "del_loc_rng_seed": 10,
        "del_len_rng_seed": 1
      },
      {
        "chromosome": [1],
        "model": "insert",
        "start_ins_frac": 0.7,
        "stop_ins_frac":  0.9,
        "phet": 0,
        "p_ins": 0.01,
        "lam_ins": 10,
        "ins_loc_rng_seed": 0,
        "ins_len_rng_seed": 1,
        "base_sel_rng_seed": 2
      }
    ]
  }

The file contains a list of parameter dictionaries, one for each instance of the model we wish to use. Models can be
repeated. To see the parameter list for each model use the explain function (or read the docs)

Seven Bridges Genomics
Current contact: kaushik.ghose@sbgenomics.com
"""
__version__ = '0.5.0'

from copy import copy
import numpy
import importlib
import json
import docopt
from variation import HOMOZYGOUS, HET1, HET2, vcopy
import variation
from fasta2wg import load_reference
import logging
logger = logging.getLogger(__name__)


def merge_variants_with_chromosome(c1, dnv):
  """
  Given an exiting chromosome (in variant format) merge new variants into it in zipper fashion
  Args:
    c1 (variant list)  - The original chromosome
    dnv (variant list) - The proposed variants
  Returns:
    c2 (variant list)  - The resultant chromosome with variant collisions arbitrated
  """
  c1_iter, dnv_iter = c1.__iter__(), dnv.__iter__()
  c2 = []
  append = c2.append
  # Try the zipper
  existing, denovo = next(c1_iter, None), next(dnv_iter, None)
  while existing is not None and denovo is not None:

    if existing.stop < denovo.POS - 1:  # Clearly non-interfering
      append(vcopy(existing))
      existing = next(c1_iter, None)
      continue

    if denovo.stop < existing.POS - 1:  # Clearly non-interfering
      append(vcopy(denovo))
      denovo = next(dnv_iter, None)
      continue

    # Potentially colliding, need to check zygosity
    if existing.het == denovo.het or \
       existing.het == HOMOZYGOUS or denovo.het == HOMOZYGOUS:  # This will collide, resolve in favor of existing
      append(vcopy(existing))
      existing, denovo = next(c1_iter, None), next(dnv_iter, None)  # Conflict resolved, both need to move on
    else:  # Won't collide, simply need to resolve who comes first
      if existing.POS <= denovo.POS:
        append(vcopy(existing))
        existing = next(c1_iter, None)
      else:
        append(vcopy(denovo))
        denovo = next(dnv_iter, None)

  # Now pick up any slack
  while existing is not None:
    append(vcopy(existing))
    existing = next(c1_iter, None)

  while denovo is not None:
    append(vcopy(denovo))
    denovo = next(dnv_iter, None)

  return c2


def merge_variants_with_genome(g1, variant_generator):
  """Given an original genome add any variants that come off the variant_generator"""
  g2 = copy(g1)
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
  return [{"model": importlib.import_module('mitty.Plugins.variants.' + k + '_plugin'), "params": v}
          for model_json in model_param_json
          for k, v in model_json.iteritems()]  # There really is only one key (the model name) and the value is the
                                               # parameter list


def apply_master_seed(models, master_seed=None):
  if master_seed is not None:
    for model in models:
      model["params"]["master_seed"] = numpy.random.RandomState(seed=int(master_seed)).randint(100000000, size=4)


def main(wg_file_name, vcf_file_name=None, param_file_name='', master_seed=None):
  """This does what the old mutate.py script used to do."""
  logger.debug('Reference file {:s}'.format(wg_file_name))
  ref = load_reference(wg_file_name)
  params = json.load(open(param_file_name, 'r'))
  models = load_variant_models(params['denovo_variant_models'])
  apply_master_seed(models, master_seed)
  g1 = create_denovo_genome(ref=ref, models=models)
  if vcf_file_name is not None:
    variation.vcf_save_gz(g1, vcf_file_name)
  return g1


if __name__ == "__main__": # pragma: no cover
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  main(cmd_args['--wg'], cmd_args['--vcf'], cmd_args['--param_file'], cmd_args['--master_seed'])
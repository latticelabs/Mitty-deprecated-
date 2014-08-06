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

import h5py
import numpy
import scipy.sparse as sparse
import importlib
import json
import docopt
from variation import HOMOZYGOUS, HET1, HET2
import variation
import logging
logger = logging.getLogger(__name__)


def initialize_mask(seq_fp, g1=None):
  """Given the sequence hdf5 file create a list of sparse arrays that represent the collision mask. This is a diplod
  mask assuming a haploid reference. Note that, for convenience we actually use 1 indexing. So our array is + 1 of the
  sequence length and element 0 is unused
  g1 is optional and indicates existing variants which should be placed on the mask. No checking is done to see if these
  variants collide (that's the caller's job)
  """
  mask = {int(c): sparse.lil_matrix((2, seq_fp['sequence/{:s}/1'.format(str(c))].size + 1), dtype='i1') for c in seq_fp['sequence']}
  if g1 is not None:
    fill_mask(mask, g1)
  return mask


def fill_mask(mask, g1):
  """Mark out the variants indicated by g1 into the mask. No checking is done to see if these variants collide
  (that's the caller's job). Modifies mask in place"""
  for chrom, variants in g1.iteritems():
    this_mask = mask[chrom]
    for v in variants:
      x0 = v.POS
      x1 = v.stop
      if v.het == HOMOZYGOUS or v.het == HET1:
        this_mask[0, x0:x1] = 1
      if v.het == HOMOZYGOUS or v.het == HET2:
        this_mask[1, x0:x1] = 1


def arbitrate_variant_collisions(g1, mask):
  """Given mask and variant data, arbitrate collisions based on the footprint. The footprint is the POS and stop
  attributes of the variant. Future expansion will use a priority attribute to arbitrate - currently variants that
  are placed earlier get absolute priority and perhaps link variants entries for complex structural variants.
  Mask is modified in place. The order of the variants determines their priority - earlier variants get preference
  in case of collisions
  """
  g2 = {}
  for chrom, variants in g1.iteritems():
    this_mask = mask[chrom]
    this_g = []
    for v in variants:
      x0 = v.POS
      x1 = v.stop
      copy_1 = this_mask[0, x0 - 1:x1 + 1].getnnz()
      copy_2 = this_mask[1, x0 - 1:x1 + 1].getnnz()  # We use a 1 base buffer around variants

      # If there are collisions, simply skip this variant
      if copy_1 and v.het != HET2:
        logger.debug('Collision at {:d}:1 {:d} - {:d}'.format(x0, x0, x1))
        continue
      if copy_2 and v.het != HET1:
        logger.debug('Collision at {:d}:2 {:d} - {:d}'.format(x0, x0, x1))
        continue

      # No collisions, valid variant, add to mask
      this_g += [v]
      if v.het != HET2:  # Either HOMOZYGOUS or HET1
        this_mask[0, x0:x1] = 1
      if v.het != HET1:  # Either HOMOZYGOUS or HET2
        this_mask[1, x0:x1] = 1
    g2[chrom] = this_g
  return g2


def add_variants_to_genome(g1, mask, variant_generator):
  """Given an original genome and its mask add any variants that come off the variant_generator"""
  for delta_g in variant_generator:
    filtered_delta_g = arbitrate_variant_collisions(delta_g, mask)
    for chrom, variants in filtered_delta_g.iteritems():
      try:
        g1[chrom].extend(variants)
      except KeyError:  # This is the first time we are seeing variants on this chromosome
        g1[chrom] = variants


def load_variant_models(param_json):
  """Given a list of models and parameters load the models in.
  {
    "variant_models": [
      {
        "model_name": {
            "param1": 1,
            "param2": "Hi"
        }
      },
      {
        "model_name": {
            "param1": 2,
            "param2": "Bye"
        }
      },
    ]
  }
  """
  return [{"model": importlib.import_module('mitty.Plugins.variants.' + k + '_plugin'), "params": v}
          for model_json in param_json['variant_models']
          for k, v in model_json.iteritems()]  # There really is only one key (the model name) and the value is the
                                               # parameter list


def create_denovo_genome(seq_fp, variant_generator_list):
  """Package functions together into a neat package."""
  mask = initialize_mask(seq_fp)
  g1 = {}
  for vg in variant_generator_list:
    add_variants_to_genome(g1, mask, vg)
  return g1


def main(wg_file_name, vcf_file_name, param_file_name, master_seed=None):
  """This does what the old mutate.py script used to do."""
  logger.debug('Reference file {:s}'.format(wg_file_name))
  with h5py.File(wg_file_name, 'r') as ref_fp:
    params = json.load(open(param_file_name, 'r'))
    models = load_variant_models(params)
    if master_seed is not None:
      for model in models:
        model["params"]["master_seed"] = numpy.random.RandomState(seed=int(master_seed))
    variant_generator_list = [model["model"].variant_generator(ref_fp, **model["params"]) for model in models]
    variation.vcf_save_gz(create_denovo_genome(ref_fp, variant_generator_list), vcf_file_name)


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  main(cmd_args['--wg'], cmd_args['--vcf'], cmd_args['--param_file'], cmd_args['--master_seed'])
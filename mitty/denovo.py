"""This module contains the machinery to use variation plugins to generate de novo variants to a given reference

This will gradually replace mutate, since it is effectively a superset of mutate. Mutate has a cooler name tho - maybe
mutate will be the command line script that wraps population functions ...



This module contains functions that generate variants with reference to a reference genome. The module can be called
as a script as well. This is useful for creating test data for MGR algorithms/data formats. The script will output
VCF file(s).

Commandline::

  Usage:
    mutate --wg=WG  --vcf=VCF  --paramfile=PFILE  [--master_seed=SEED] [-v]

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

import subprocess
import pysam
import tempfile
import h5py
import numpy
import scipy.sparse as sparse
import os
import imp
import json
import docopt
import datetime
import logging
from variation import HOMOZYGOUS, HET1, HET2
logger = logging.getLogger(__name__)


def load_models(json_list_models):
  """Given a list of models and parameters load the models in."""
  plugin_dir = os.path.join(os.path.dirname(__file__), 'Plugins', 'variants')  # Thanks Nebojsa Tijanic!
  models = []
  for model_json in json_list_models:
    fp, pathname, description = imp.find_module(model_json['model'] + '_plugin', [plugin_dir])
    models.append(imp.load_module(model_json['model'], fp, pathname, description))
  return models


def initialize_mask(seq_fp):
  """Given the sequence hdf5 file create a list of sparse arrays that represent the collision mask. This is a diplod
  mask assuming a haploid reference. Note that, for convenience we actually use 1 indexing. So our array is + 1 of the
  sequence length and element 0 is unused"""
  return {c: sparse.lil_matrix((2, seq_fp['sequence/{:s}/1'.format(str(c))].size + 1), dtype='i1') for c in seq_fp['sequence']}


def arbitrate_variant_collisions(mask, g1):
  """Given mask and variant data, arbitrate collisions based on the footprint. The footprint is the POS and stop
  attributes of the variant. Future expansion will use a priority attribute to arbitrate - currently variants that
  are placed earlier get absolute priority.
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


def filter_variants(mask, variants):
  """Given mask and variant data, arbitrate collisions based on the footprint. The footprint is the POS and

  footprint   -  [ [(het, chrom, pos_st, pos_nd) ...] ... ]


  >>> f = h5py.File("Examples/Data/chimera.h5","r")
  >>> mask = initialize_mask(f)
  >>> footprint = [[(1, 0, 100, 101)]]  # This should pass as a variant
  >>> filter_variants(mask, footprint)
  [0]
  >>> footprint = [[(1, 0, 70, 121)]]  # This should not pass as a variant
  >>> filter_variants(mask, footprint)
  []
  >>> footprint = [[(2, 0, 70, 121)]]  # This should pass as a variant - collision is on other chromosome copy
  >>> filter_variants(mask, footprint)
  [0]
  >>> footprint = [[(1, 0, 70, 80), (1, 0, 90, 110)]]  # This should not pass as a variant
  >>> filter_variants(mask, footprint)
  []
  >>> footprint = [[(1, 0, 70, 80)], [(1, 0, 90, 110)]]  # Only the first should pass as a variant
  >>> filter_variants(mask, footprint)
  [0]
  >>> footprint = [[(1, 0, 120, 130)], [(2, 0, 140, 150)]]  # All should pass as variants
  >>> filter_variants(mask, footprint)
  [0, 1]
  >>> footprint = [[(3, 0, 90, 110)], [(3, 0, 140, 150)]]  # None should pass as variants
  >>> filter_variants(mask, footprint)
  []
  """
  idx = []
  for n, feet in enumerate(footprint):
    collision = False
    for foot in feet:
      het = foot[0]
      chrom = foot[1]
      copy_1 = mask[chrom][0, foot[2] - 1:foot[3] + 1].getnnz()  # footprint is internal and is 0 indexed
      copy_2 = mask[chrom][1, foot[2] - 1:foot[3] + 1].getnnz()  # We use a 1 base buffer around variants
      # If there are collisions, simply skip this variant
      if copy_1:
        if het == 1 or het == 3:
          collision = True
          logger.debug('Collision at {:d}:1 {:d} - {:d}'.format(foot[1] + 1, foot[2], foot[3]))
          break
      if copy_2:
        if het == 2 or het == 3:
          collision = True
          logger.debug('Collision at {:d}:2 {:d} - {:d}'.format(foot[1] + 1, foot[2], foot[3]))
          break

    if not collision:
      # No collisions, valid variant, add to mask
      idx.append(n)
      if het == 1 or het == 3:
        mask[chrom][0, foot[2]:foot[3]] = 1
      if het == 2 or het == 3:
        mask[chrom][1, foot[2]:foot[3]] = 1
  return idx


def write_vcf_header(file_handle, sim_date, argv, reference_filename):
  """Given a file handle, write out a suitable header to start the VCF file
  Inputs:
    file_handle       - an open file handle
    sim_date          - date in string format. Get's dumped into 'fileDate'
    argv              - get's dumped into the 'source' string of the vcf file
    reference_filename  - name of the reference, dumped into the 'reference' string

  Notes:
  1. 'filedate' is the date of the simulation
  2. 'source' contains the version of the mutate program used to generate the data
  3. 'reference' contains the name of the file entered as the reference file

  """
  file_handle.write(
    """##fileformat=VCFv4.1
##fileDate={:s}
##source=mutate.py {:s} ({:s})
##reference={:s}
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n""".
    format(sim_date, __version__, argv, reference_filename)
  )


def write_vcf_mutations(fp, variants):
  """Given a mutator format dictionary write the mutations in VCF format into the file
  Inputs:
    file_handle   - handle of an opened text file. The output will be appended to this file.
    chrom         - chromosome number
    variants      - list of tuples (POS, REF, ALT) : standard format as returned by the variant plugins
  """
  for var in variants:
    #       CHROM    POS   ID   REF   ALT   QUAL FILTER INFO FORMAT tsample
    fp.write("{:d}\t{:d}\t{:s}\t{:s}\t{:s}\t{:d}\t{:s}\t{:s}\t{:s}\t{:s}\n".format(*var))


def create_mutations(wg_file_name, params, master_seed=None):
  """Same as create_and_save_mutations, except that mutations are held in memory, making this method faster e.g. for
  population operations"""
  logger.debug('Reference file {:s}'.format(wg_file_name))

  with h5py.File(wg_file_name, 'r') as ref_fp, open(vcf_file_name, 'w') as vcf_fp:
    model_list = load_models(params['variant_models'])
    if master_seed is None:
      ms_rng = None
    else:
      ms_rng = numpy.random.RandomState(seed=int(master_seed))
    mask = initialize_mask(ref_fp)

    write_vcf_header(vcf_fp, datetime.datetime.now().isoformat(), docopt.sys.argv.__str__(),
                     os.path.basename(wg_file_name))
    for model, model_params in zip(model_list, params['variant_models']):
      logger.debug('Running {:s}'.format(model_params['model']))
      if ms_rng:
        model_params['master_seed'] = ms_rng.randint(100000000)
        logger.debug('Sending master seed of {:d}'.format(model_params['master_seed']))

      description, footprint, vcf_line = model.variants(ref_fp, **model_params)
      logger.debug('{:d} {:s} generated'.format(len(vcf_line), model_params['model']))
      idx = filter_variants(mask, footprint)
      logger.debug('{:d} {:s} placed without collision'.format(len(idx), model_params['model']))
      write_vcf_mutations(vcf_fp, [vl for n in idx for vl in vcf_line[n]])


def create_and_save_mutations(vcf_file_name, wg_file_name, params, master_seed=None):
  """This routine saves variants in batches and uses less memory. For faster, in memory version of this routine
  look at create_mutations"""
  logger.debug('Reference file {:s}\nVCF file {:s}'.format(wg_file_name, vcf_file_name))

  with h5py.File(wg_file_name, 'r') as ref_fp, open(vcf_file_name, 'w') as vcf_fp:
    model_list = load_models(params['variant_models'])
    if master_seed is None:
      ms_rng = None
    else:
      ms_rng = numpy.random.RandomState(seed=int(master_seed))
    mask = initialize_mask(ref_fp)

    write_vcf_header(vcf_fp, datetime.datetime.now().isoformat(), docopt.sys.argv.__str__(),
                     os.path.basename(wg_file_name))
    for model, model_params in zip(model_list, params['variant_models']):
      logger.debug('Running {:s}'.format(model_params['model']))
      if ms_rng:
        model_params['master_seed'] = ms_rng.randint(100000000)
        logger.debug('Sending master seed of {:d}'.format(model_params['master_seed']))

      description, footprint, vcf_line = model.variants(ref_fp, **model_params)
      logger.debug('{:d} {:s} generated'.format(len(vcf_line), model_params['model']))
      idx = filter_variants(mask, footprint)
      logger.debug('{:d} {:s} placed without collision'.format(len(idx), model_params['model']))
      write_vcf_mutations(vcf_fp, [vl for n in idx for vl in vcf_line[n]])


def sort_vcf(in_vcf_name, out_vcf_name):
  #vcf-sort the.vcf > sorted.vcf
  logger.debug('Sorting {:s}'.format(in_vcf_name))
  with open(out_vcf_name, 'w') as fp:
    subprocess.call(['vcf-sort', in_vcf_name], stdout=fp)


def compress_and_index_vcf(in_vcf_name, out_vcf_name):
  #bgzip -c sorted.vcf > sorted.vcf.gz
  #tabix sorted.vcf.gz
  logger.debug('Compressing and indexing {:s}'.format(in_vcf_name))
  pysam.tabix_compress(in_vcf_name, out_vcf_name, force=True)
  pysam.tabix_index(out_vcf_name, force=True, preset='vcf')


def main(args):
  """variants take the form of ."""
  vcf_gz_name = args['--vcf']
  vcf_name, ext = os.path.splitext(vcf_gz_name)
  if ext != '.gz':
    vcf_name += '_srt.vcf'

  temp_vcf_fp, temp_vcf_name = tempfile.mkstemp(suffix='.vcf')
  params = json.load(open(args['--paramfile'], 'r'))
  create_and_save_mutations(temp_vcf_name, args['--wg'], params, args['--master_seed'])
  sort_vcf(temp_vcf_name, vcf_name)
  compress_and_index_vcf(vcf_name, vcf_gz_name)
  os.remove(temp_vcf_name)


if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  else:
    cmd_args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  main(cmd_args)
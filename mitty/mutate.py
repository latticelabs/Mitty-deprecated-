"""This module contains functions that generate variants with reference to a reference genome. The module can be called
as a script as well. This is useful for creating test data for MGR algorithms/data formats. The script will output
VCF file(s).

Commandline::

  Usage:
    mutate --wg=WG  --vcf=VCF  --paramfile=PFILE  [--master_seed=SEED] [-v]

  Options:
    --wg=WG                 Whole genome reference file
    --vcf=VCF               Output VCF file
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
__version__ = '0.4.0'

import h5py
import numpy
import scipy.sparse as sparse
import os
import imp
import json
import docopt
import datetime
import logging

logger = logging.getLogger(__name__)


def load_models(json_list_models):
  plugin_dir = os.path.join(os.path.dirname(__file__), 'Plugins', 'Mutation')  # Thanks Nebojsa Tijanic!
  models = []
  for model_json in json_list_models:
    #model_fname = os.path.join(plugin_dir, model_json['model'] + '_plugin.py')
    #models.append(imp.load_source(model_fname, model_fname, open(model_fname, 'r')))
    fp, pathname, description = imp.find_module(model_json['model'] + '_plugin', [plugin_dir])
    models.append(imp.load_module(model_json['model'], fp, pathname, description))
  return models


def initialize_mask(seq_fp):
  """Given the sequence hdf5 file create a list of sparse arrays that represent the collision mask. This is a diplod
  mask assuming a haploid reference """
  return [sparse.lil_matrix((2, seq_fp['sequence/{:s}/1'.format(str(c))].size), dtype='i1') for c in seq_fp['sequence']]


def filter_variants(mask, footprint):
  """Given mask and variant data, arbitrate collisions.

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


def main(args):
  """variants take the form of ."""
  wg_file_name = args['--wg']
  vcf_file_name = args['--vcf']

  with h5py.File(wg_file_name, 'r') as ref_fp, open(vcf_file_name, 'w') as vcf_fp:

    params = json.load(open(args['--paramfile'], 'r'))
    model_list = load_models(params['variant_models'])
    ms_rng = numpy.random.RandomState(seed=int(args['--master_seed'])) if args['--master_seed'] is not None else None
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

  #For further use, this vcf file usually needs to be sorted and then compressed and indexed
  #vcf-sort the.vcf > sorted.vcf
  #bgzip -c sorted.vcf > sorted.vcf.gz
  #tabix sorted.vcf.gz

  # logger.debug('Sorting VCF file')
  # print v_file_name
  # ps = subprocess.Popen(('vcf-sort', v_file_name), stdout=subprocess.PIPE)
  # output = subprocess.check_output(('temp.vcf'), stdin=ps.stdout)
  # ps.wait()
  #
  # #subprocess.call(['vcf-sort', v_file_name, '>', 'temp.vcf'], shell=True)
  # logger.debug('Compressing and indexing VCF file')
  # pysam.tabix_compress('temp.vcf', v_file_name + '.gz', force=True)
  # pysam.tabix_index(v_file_name + '.gz', force=True, preset='vcf')

if __name__ == "__main__":
  if len(docopt.sys.argv) < 2:  # Print help message if no options are passed
    docopt.docopt(__doc__, ['-h'])
  elif docopt.sys.argv[1] == 'explain':
    print __explain__
    exit(0)
  else:
    cmd_args = docopt.docopt(__doc__, version=__version__)

  level = logging.DEBUG if cmd_args['-v'] else logging.WARNING
  logging.basicConfig(level=level)

  main(cmd_args)
import tempfile
import json

import mitty.lib.io as mio
import mitty.lib.variants as vr
import mitty.genomes as genomes
import mitty.plugins.variants.snp_plugin as snp
import mitty.plugins.site_frequency.double_exp as double_exp
import mitty.plugins.population.standard as standard
from mitty.tests import *  # To get definitions from the setup script


def run_simulations_test():
  """Core genome simulation loop"""
  _, temp_name = tempfile.mkstemp(suffix='.hdf5')

  pop_db_name = temp_name
  ref = mio.Fasta(multi_dir=example_data_dir)
  sfs_model = double_exp.Model()
  variant_models = [snp.Model(p=0.001), snp.Model(p=0.005)]
  chromosomes = [1, 2]
  sample_size = 20
  master_seed = 2
  pop_model = standard.Model(sample_size=sample_size)

  genomes.run_simulations(pop_db_name, ref=ref, sfs_model=sfs_model, variant_models=variant_models,
                          population_model=pop_model, chromosomes=chromosomes, master_seed=master_seed)

  # The parts are tested separately, we just want a superficial test to see if we wrote the correct number of things
  pop = vr.Population(fname=pop_db_name)
  ml = pop.get_master_list(chrom=1)
  index_lists = [[pop.get_sample_chromosome(chrom=ch, sample_name='g{:d}_s{:d}'.format(0, n)) for ch in [1, 2]] for n in range(sample_size)]
  null_ml = pop.get_master_list(chrom=3)
  null_index_list = pop.get_sample_chromosome(chrom=3, sample_name='g0_s0')
  os.remove(temp_name)

  for n in range(sample_size):
    for ch in [0, 1]:
      assert index_lists[n][ch].shape[0] > 0
  assert len(null_ml) == 0
  assert null_index_list.shape[0] == 0


def integration_test():
  """'genomes' command line program"""
  _, param_file = tempfile.mkstemp(suffix='.json')
  _, db_file = tempfile.mkstemp(suffix='.hdf5')
  test_params = {
    "files": {
      "reference_dir": example_data_dir,
      "dbfile": db_file
    },
    "rng": {
      "master_seed": 1
    },
    "sample_size": 10,
    "site_model": {
        "double_exp": {
          "k1": 0.1,
          "k2": 2.0,
          "p0": 0.001,
          "p1": 0.2,
          "bin_cnt": 30
        }
    },
    "chromosomes": [1, 2],
    "variant_models": [
      {
        "snp": {
          "p": 0.001
        }
      }
    ]
  }

  json.dump(test_params, open(param_file, 'w'))
  genomes.generate({'<pfile>': param_file, '-p': False})

  pop = vr.Population(fname=db_file)
  ml = pop.get_master_list(chrom=1)

  os.remove(param_file)
  os.remove(db_file)

  assert len(ml) > 0
  assert len(pop.get_sample_names()) == 10
  assert 'g0_s6' in pop.get_sample_names()
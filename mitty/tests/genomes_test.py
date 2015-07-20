import tempfile
import json

from click.testing import CliRunner

import mitty.lib.variants as vr
import mitty.genomes as genomes
from mitty.tests import *  # To get definitions from the setup script


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

  runner = CliRunner()
  result = runner.invoke(genomes.cli, ['generate', param_file])
  assert result.exit_code == 0, result
  assert os.path.exists(db_file)

  pop = vr.Population(fname=db_file)
  ml = pop.get_master_list(chrom=1)

  os.remove(param_file)
  os.remove(db_file)

  assert len(ml) > 0
  assert len(pop.get_sample_names()) == 10
  assert 'g0_s6' in pop.get_sample_names()
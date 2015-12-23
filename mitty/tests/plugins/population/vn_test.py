"""Run genomes generate and then analyze the resultant genomes file to make sure it satisfies the conditions we
set it"""
import tempfile
import os
import json

from click.testing import CliRunner

import mitty.lib.variants as vr
import mitty.genomes as genomes
import mitty.tests


def vn_test():
  """Test 'vn' population model"""
  _, param_file = tempfile.mkstemp(suffix='.json')
  _, db_file = tempfile.mkstemp(suffix='.hdf5')
  test_params = {
    "files": {
      "reference_file": mitty.tests.test_fasta_genome_file,
      "dbfile": db_file
    },
    "rng": {
      "master_seed": 12345
    },
    "population_model": {
      "vn": {
        "p_vx": 0.2,
        "p_vn": [0.1, 0.5, 0.9]
      }
    },
    "chromosomes": [1],
    "variant_models": [
      {
        "snp": {
          "p": 0.01
        }
      }
    ]
  }
  json.dump(test_params, open(param_file, 'w'))

  runner = CliRunner()
  result = runner.invoke(genomes.cli, ['generate', param_file])
  assert result.exit_code == 0, result
  assert os.path.exists(db_file)

  pop = vr.Population(fname=db_file, mode='r', in_memory=False)
  # ml = pop.get_variant_master_list(chrom=4)

  chrom_idx_vx = pop.get_sample_variant_index_for_chromosome(1, 'vx')
  idx_vx = [set([i[0] for i in chrom_idx_vx if i[1] != c]) for c in [1, 0]]

  idx_vn = []
  for v in ['v0', 'v1', 'v2']:
    chrom_idx = pop.get_sample_variant_index_for_chromosome(1, v)
    idx_vn.append([set([i[0] for i in chrom_idx if i[1] != c]) for c in [1, 0]])

  for n in range(len(idx_vn) - 1):
    for cpy in [0, 1]:
      assert len(idx_vx[cpy]) > 0, idx_vx[cpy]
      assert abs(len(idx_vn[n][cpy].intersection(idx_vn[n + 1][cpy])) - len(idx_vn[n][cpy])) < 0.05 * len(idx_vn[n][cpy])
      assert abs(len(idx_vx[cpy].intersection(idx_vn[n][cpy])) - len(idx_vx[cpy].intersection(idx_vn[n + 1][cpy]))) < 0.05 * len(idx_vx[cpy])

  os.remove(param_file)
  os.remove(db_file)

# Need cleanup code to work even if test fails ...
# nosetests mitty.tests.plugins.population
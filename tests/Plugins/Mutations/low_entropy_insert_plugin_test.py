import json
import mitty.mutate
from ... import *


def integration_test():
  """Testing if low_entropy_insert works with mutate"""
  vcf_name = os.path.join(data_dir , 'low_entropy_test.vcf.gz')
  param_file = os.path.join(data_dir, 'low_entropy_test.json')

  params = {
    "variant_models": [
      {
        "chromosome": [2],
        "model": "low_entropy_insert",
        "phet": 0.9,
        "p": 0.01,
        "ins_len_lo": 10,
        "ins_len_hi": 20,
        "sub_seq_len": 5,
        "copy_from_neighborhood": False,
        "master_seed": 0,
        "ins_loc_rng_seed": 1,
        "ins_len_rng_seed": 2,
        "base_sel_rng_seed": 3,
        "het_rng_seed": 4,
        "copy_rng_seed": 5
      },
      {
        "chromosome": [2],
        "model": "low_entropy_insert",
        "phet": 0.9,
        "p": 0.01,
        "ins_len_lo": 10,
        "ins_len_hi": 20,
        "sub_seq_len": 10,
        "copy_from_neighborhood": True,
        "master_seed": 0,
        "ins_loc_rng_seed": 1,
        "ins_len_rng_seed": 2,
        "base_sel_rng_seed": 3,
        "het_rng_seed": 4,
        "copy_rng_seed": 5
      }

    ]
  }
  json.dump(params, open(param_file, 'w'), indent=2)

  args = {
    '--wg': wg_name,
    '--vcf': vcf_name,
    '--paramfile': param_file,
    '--master_seed': 0,
    '-v': True
  }
  mitty.mutate.main(args)

  assert os.path.exists(vcf_name)  # A very simple test

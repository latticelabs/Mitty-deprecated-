"""The snp_plugin does not run by itself (it merely runs tests), but this wrapper is used to generate parts of the
.json file that will be used by mutate_wrapper to create a full .json file that can be fed to mutate.py

The .json fragment should look like
          ______________ model id
        /
      "snp": {
          "model": "snp",
          "start_snps_frac": 0.1,
          "stop_snps_frac":  0.3,
          "phet": 0.5,
          "p": 0.01,
          "het_rng_seed":3,
          "strand_rng_seed": 4,
          "poisson_rng_seed": 1,
          "base_sub_rng_seed": 2
      }

mutate_wrapper will place this snippet under "mutations"
"""
import os
import json
from sbgsdk import define
from nose.tools import assert_equals

class SNP(define.Wrapper):
  class Inputs(define.Inputs):
    pass
    # ref = define.input(name='Reference', description='The reference file')

  class Outputs(define.Outputs):
    json_fragment = define.output(name="SNPs", description="SNP definitions for mutate")

  class Params(define.Params):
    model_id = define.string(required=True, description='A unique name for this instance of the SNP generator',
                             category='General')
    # "snp" in the example above. Needs to be unique
    start_snps_frac = define.real(default=0, min=0, max=1, category='Model params',
                                  description='start generating snps from here (0.0, 1.0)')
    stop_snps_frac = define.real(default=1, min=0, max=1, category='Model params',
                                 description='stop generating snps after this (0.0, 1.0)')
    phet = define.real(default=0.01, min=0, max=1, category='Model params',
                       description='probability of having heterozygous mutation')
    p = define.real(default=0.01, min=0, max=1, category='Model params',
                    description='probability of SNPs')
    het_rng_seed = define.integer(default=1, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for random number generator used to decide if genotype is heterozygous or not')
    strand_rng_seed = define.integer(default=1, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for random number generator used to decide which strand the SNP will be on')
    poisson_rng_seed = define.integer(default=1, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for SNP locator random number generator')
    base_sub_rng_seed = define.integer(default=1, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for random number generator used to select ALT bases')

  def write_to_json(self, fname):
    with open(fname, 'w') as f:
      params = self.params.__json__()
      params.pop('model_id')
      json.dump({self.params.model_id: dict(model='snp', **params)}, f)

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    self.outputs.json_fragment = \
      os.path.join(output_dir, '{:s}_snp_plugin_params.json'.format(self.params.model_id))
    # By adding the model_id bit to the name we ensure uniqueness
    self.write_to_json(self.outputs.json_fragment)
    self.outputs.json_fragment.meta = self.outputs.json_fragment.make_metadata(file_type='json')


def test():
  params = {
            "model_id": "snp_test",
            "start_snps_frac": 0.1,
            "stop_snps_frac": 0.3,
            "phet": 0,
            "p": 0.01,
            "het_rng_seed": 3,
            "strand_rng_seed": 4,
            "poisson_rng_seed": 1,
            "base_sub_rng_seed": 2
        }
  inputs = {}
  wrp = SNP(inputs, params)
  outputs = wrp.test()
  expected = {"snp_test": {
            "model": "snp",
            "start_snps_frac": 0.1,
            "stop_snps_frac":  0.3,
            "phet": 0,
            "p": 0.01,
            "het_rng_seed":3,
            "strand_rng_seed": 4,
            "poisson_rng_seed": 1,
            "base_sub_rng_seed": 2
        }}
  with open(outputs.json_fragment) as fp:
    assert_equals(json.load(fp), expected)

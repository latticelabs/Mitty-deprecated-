"""The snp_plugin does not run by itself (it merely runs tests), but this wrapper is used to generate parts of the
.json file that will be used by mutate_wrapper to create a full .json file that can be fed to mutate.py

The .json fragment should look like

    "snp": {
        "model": "snp",
        "start_snps_frac": 0.1,
        "stop_snps_frac":  0.3,
        "phet": 0,
        "p": 0.01,
        "poisson_rng_seed": 1,
        "base_sub_rng_seed": 2
    },

mutate_wrapper will place this snippet under "mutations"
"""
from sbgsdk import define, Process, require
import json
from nose.tools import assert_equals


class SNP(define.Wrapper):
  class Inputs(define.Inputs):
    pass
    # ref = define.input(name='Reference', description='The reference file')

  class Outputs(define.Outputs):
    json_fragment = define.output()

  class Params(define.Params):
    model_id = define.string(required=True)   # "snp" in the example above.
                                              # Needs to be unique8, currently not used by mutate
    start_snps_frac = define.real(default=0, min=0, max=1, category='SNP model')
    stop_snps_frac = define.real(default=1, min=0, max=1, category='SNP model')
    phet = define.real(default=0.01, min=0, max=1, category='SNP model')
    p = define.real(default=0.01, min=0, max=1, category='SNP model')
    poisson_rng_seed = define.integer(default=1, min=0, category='SNP model')
    base_sub_rng_seed = define.integer(default=1, min=0, category='SNP model')

    # def make_json(self):
    #   params = self.__json__()
    #   params.pop('model_id')
    #   return {self.model_id: dict(model='snp', **params)}

  def execute(self):
    self.outputs.json_fragment = 'snp_plugin_params.json'
    with open(self.outputs.json_fragment, 'w') as f:
      params = self.params.__json__()
      params.pop('model_id')
      json.dump({self.params.model_id: dict(model='snp', **params)}, f)
      #json.dump(self.params.make_json(), f)


def test():
  params = {
            "model_id": "snp_test",
            "start_snps_frac": 0.1,
            "stop_snps_frac":  0.3,
            "phet": 0,
            "p": 0.01,
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
            "poisson_rng_seed": 1,
            "base_sub_rng_seed": 2
        }}
  with open(outputs.json_fragment) as fp:
    assert_equals(json.load(fp), expected)
  #assert os.path.getsize('snp_plugin_params.json')
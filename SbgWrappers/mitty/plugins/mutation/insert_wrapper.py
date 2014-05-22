"""This wrapper only translates input parameters into a json file that is an input to the mutate wrapper. The mutate
wrapper combines json fragments from all the plugins into a complete json parameters file that is then fed to mutate.py

The .json fragment should look like

          ______________ model id
        /
    "myinsert": {
        "model": "insert",
        "start_ins_frac": 0.7,
        "stop_ins_frac":  0.9,
        "phet": 0,
        "p_ins": 0.01,
        "lam_ins": 10,
        "het_rng_seed": 3,
        "strand_rng_seed": 4,
        "ins_loc_rng_seed": 0,
        "ins_len_rng_seed": 1,
        "base_sel_rng_seed": 2
    }

mutate_wrapper will place this snippet under "mutations"
"""
import os
import json
from sbgsdk import define
from nose.tools import assert_equals


class Insert(define.Wrapper):
  class Inputs(define.Inputs):
    pass

  class Outputs(define.Outputs):
    json_fragment = define.output(name="Inserts", description="Insertion definitions for mutate")

  class Params(define.Params):
    model_id = define.string(required=True,
                             description='A unique name for this instance of the insert generator', category='General')
    # "myinsert" in the example above. Needs to be unique TODO: How to make interlock on platform
    start_ins_frac = define.real(default=0, min=0, max=1, category='Model params',
                                  description='start generating inserts from here (0.0, 1.0)')
    stop_ins_frac = define.real(default=1, min=0, max=1, category='Model params',
                                 description='stop generating inserts after this (0.0, 1.0)')
    phet = define.real(default=0.01, min=0, max=1, category='Model params',
                       description='probability of having heterozygous mutation')
    p_ins = define.real(default=0.01, min=0, max=1, category='Model params',
                        description='probability of an insert')
    lam_ins = define.integer(default=10, min=1, category='Model params',
                             description='Mean length of insertion (Poisson distributed)')
    het_rng_seed = define.integer(default=1, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for random number generator used to decide if genotype is heterozygous or not')
    strand_rng_seed = define.integer(default=2, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for random number generator used to decide which strand the variant will be on')
    ins_loc_rng_seed = define.integer(default=3, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for variant location random number generator')
    ins_len_rng_seed = define.integer(default=4, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for insert length random number generator')
    base_sel_rng_seed = define.integer(default=5, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for random number generator used to select ALT bases')

  def write_to_json(self, fname):
    with open(fname, 'w') as f:
      params = self.params.__json__()
      params.pop('model_id')
      json.dump({self.params.model_id: dict(model='insert', **params)}, f)

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    self.outputs.json_fragment = \
      os.path.join(output_dir, '{:s}_insert_plugin_params.json'.format(self.params.model_id))
    # By adding the model_id bit to the name we ensure uniqueness
    self.write_to_json(self.outputs.json_fragment)
    self.outputs.json_fragment.meta = self.outputs.json_fragment.make_metadata(file_type='json')


def test():
  params = {
      "model_id": "insert_test",
      "start_ins_frac": 0.7,
      "stop_ins_frac":  0.9,
      "phet": 0,
      "p_ins": 0.01,
      "lam_ins": 10,
      "het_rng_seed": 3,
      "strand_rng_seed": 4,
      "ins_loc_rng_seed": 0,
      "ins_len_rng_seed": 1,
      "base_sel_rng_seed": 2
  }
  inputs = {}
  wrp = Insert(inputs, params)
  outputs = wrp.test()
  expected = {'insert_test': {
        "model": "insert",
        "start_ins_frac": 0.7,
        "stop_ins_frac":  0.9,
        "phet": 0,
        "p_ins": 0.01,
        "lam_ins": 10,
        "het_rng_seed": 3,
        "strand_rng_seed": 4,
        "ins_loc_rng_seed": 0,
        "ins_len_rng_seed": 1,
        "base_sel_rng_seed": 2
    }
  }
  with open(outputs.json_fragment) as fp:
    assert_equals(json.load(fp), expected)

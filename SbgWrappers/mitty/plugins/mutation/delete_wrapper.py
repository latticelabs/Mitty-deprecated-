"""This wrapper only translates input parameters into a json file that is an input to the mutate wrapper. The mutate
wrapper combines json fragments from all the plugins into a complete json parameters file that is then fed to mutate.py

The .json fragment should look like

          ______________ model id
        /
    "mydelete": {
        "model": "delete",
        "start_dels_frac": 0.7,
        "stop_dels_frac":  0.9,
        "phet": 0,
        "p_del": 0.01,
        "lam_del": 10,
        "het_rng_seed": 3,
        "strand_rng_seed": 4,
        "del_loc_rng_seed": 0,
        "del_len_rng_seed": 1
    }

mutate_wrapper will place this snippet under "mutations"
"""
import os
import json
from sbgsdk import define
from nose.tools import assert_equals


class Deletion(define.Wrapper):
  class Inputs(define.Inputs):
    pass

  class Outputs(define.Outputs):
    json_fragment = define.output(name="Deletions", description="Deletion definitions for mutate")

  class Params(define.Params):
    model_id = define.string(required=True,
                             description='A unique name for this instance of the insert generator', category='General')
    # "mydelete" in the example above. Needs to be unique TODO: How to make interlock on platform
    start_dels_frac = define.real(default=0, min=0, max=1, category='Model params',
                                  description='start generating deletions from here (0.0, 1.0)')
    stop_dels_frac = define.real(default=1, min=0, max=1, category='Model params',
                                 description='stop generating deletions after this (0.0, 1.0)')
    phet = define.real(default=0.01, min=0, max=1, category='Model params',
                       description='probability of having heterozygous mutation')
    p_del = define.real(default=0.01, min=0, max=1, category='Model params',
                        description='probability of a deletion')
    lam_del = define.integer(default=10, min=1, category='Model params',
                             description='Mean length of deletion (Poisson distributed)')
    het_rng_seed = define.integer(default=1, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for random number generator used to decide if genotype is heterozygous or not')
    strand_rng_seed = define.integer(default=2, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for random number generator used to decide which strand the data will be on')
    del_loc_rng_seed = define.integer(default=3, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for data location random number generator')
    del_len_rng_seed = define.integer(default=4, min=0, max=2**32 - 1, category='Model params: RNG',
          description='Seed for deletion length random number generator')

  def write_to_json(self, fname):
    with open(fname, 'w') as f:
      params = self.params.__json__()
      params.pop('model_id')
      json.dump({self.params.model_id: dict(model='delete', **params)}, f)

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    self.outputs.json_fragment = \
      os.path.join(output_dir, '{:s}_delete_plugin_params.json'.format(self.params.model_id))
    # By adding the model_id bit to the name we ensure uniqueness
    self.write_to_json(self.outputs.json_fragment)
    self.outputs.json_fragment.meta = self.outputs.json_fragment.make_metadata(file_type='json')


def test():
  params = {
    "model_id": "delete_test",
    "model": "delete",
    "start_dels_frac": 0.7,
    "stop_dels_frac":  0.9,
    "phet": 0,
    "p_del": 0.01,
    "lam_del": 10,
    "het_rng_seed": 3,
    "strand_rng_seed": 4,
    "del_loc_rng_seed": 0,
    "del_len_rng_seed": 1
  }
  inputs = {}
  wrp = Deletion(inputs, params)
  outputs = wrp.test()
  expected = {"delete_test": {
        "model": "delete",
        "start_dels_frac": 0.7,
        "stop_dels_frac":  0.9,
        "phet": 0,
        "p_del": 0.01,
        "lam_del": 10,
        "het_rng_seed": 3,
        "strand_rng_seed": 4,
        "del_loc_rng_seed": 0,
        "del_len_rng_seed": 1
    }
  }
  with open(outputs.json_fragment) as fp:
    assert_equals(json.load(fp), expected)

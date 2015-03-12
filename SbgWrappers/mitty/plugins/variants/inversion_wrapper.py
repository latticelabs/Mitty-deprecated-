"""
This wrapper converts the parameters into a json file that is fed into the denovo wrap.

Schema example:

    {
        "chromosome": [1],
        "model": "inversion",
        "phet": 0.5,
        "p": 0.01,
        "inv_len_lo": 10,
        "inv_len_hi": 100,
        "inv_loc_rng_seed": 1,
        "inv_len_rng_seed": 2,
        "het_rng_seed": 3,
        "copy_rng_seed": 4
    }
"""

import os
import json
from sbgsdk import define
from nose.tools import assert_equals

class Inversion(define.Wrapper):
  class Inputs(define.Inputs):
    pass

  class Outputs(define.Outputs):
    json_fragment = define.output(name="Inversion", description="Inversion definitions for denovo")

  class Params(define.Params):
  	model_id = define.string(required=True,
                             description='A unique name for this instance of the inversion generator', category='General')
  	chromosome = define.integer(default=1, list=True)
  	p = define.real(default=0.01, min=0, max=1, category='Model params', description='probability of an inversion')
  	phet = define.real(default=0.05, min=0, max=1, category='Model params',
  	                   description='probability of having heterozygous mutation')
  	inv_len_lo = define.integer(default=10, min=1, category='Model params',
  	                         description='Lower bound on inversion lengths')
  	inv_len_hi = define.integer(default=100, min=1, category='Model params',
  	                         description='Upper bound on inversion lengths')
  	inv_loc_rng_seed = define.integer(default=1, min=0, max=2**32 - 1, category='Model params: RNG',
  	      description='Seed for data location random number generator')
  	inv_len_rng_seed = define.integer(default=2, min=0, max=2**32 - 1, category='Model params: RNG',
  	      description='Seed for inversion length random number generator')
  	het_rng_seed = define.integer(default=3, min=0, max=2**32 - 1, category='Model params: RNG',
  	      description='Seed for inversion length random number generator')
  	copy_rng_seed = define.integer(default=4, min=0, max=2**32 - 1, category='Model params: RNG')

  def write_to_json(self, fname):
  	with open(fname, 'w') as f:
  	  params = self.params.__json__()
  	  params.pop('model_id')
  	  json.dump({'inversion': dict(**params)}, f)

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    self.outputs.json_fragment = \
      os.path.join(output_dir, '{:s}_inversion_plugin_params.json'.format(self.params.model_id))
    # By adding the model_id bit to the name we ensure uniqueness
    self.write_to_json(self.outputs.json_fragment)
    self.outputs.json_fragment.meta = self.outputs.json_fragment.make_metadata(file_type='json')

def test():
  params = {
      "chromosome": [1],
      "model_id": "inversion",
      "phet": 0.5,
      "p": 0.01,
      "inv_len_lo": 10,
      "inv_len_hi": 100,
      "inv_loc_rng_seed": 1,
      "inv_len_rng_seed": 2,
      "het_rng_seed": 3,
      "copy_rng_seed": 4
  }
  inputs = {}
  wrp = Inversion(inputs, params)
  outputs = wrp.test()
  expected = {"inversion": {
      "chromosome": [1],
      "phet": 0.5,
      "p": 0.01,
      "inv_len_lo": 10,
      "inv_len_hi": 100,
      "inv_loc_rng_seed": 1,
      "inv_len_rng_seed": 2,
      "het_rng_seed": 3,
      "copy_rng_seed": 4
	    }
  }
  assert_equals.__self__.maxDiff = None
  with open(outputs.json_fragment) as fp:
    assert_equals(json.load(fp), expected)
"""
This wrapper converts the parameters into a json file that is fed into the denovo wrap.

Schema example:

{
  "chromosome": [1],      # List of chromosomes to apply the variant to
  "p": 0.01,              # probability that the deletion will happen at any given base
  "phet": 0.5,            # probability that the variant will be heterozygous
  "del_len_lo": 100,      # lower bound on deletion lengths
  "del_len_hi": 10000,    # upper bound on deletion lengths
  "del_loc_rng_seed": 1,  #
  "del_len_rng_seed": 2,
  "het_rng_seed": 3,
  "copy_rng_seed": 4
}
"""

import os
import json
from sbgsdk import define
from nose.tools import assert_equals

class Deletion(define.Wrapper):
  class Inputs(define.Inputs):
    pass

  class Outputs(define.Outputs):
    json_fragment = define.output(name="Deletions", description="Deletion definitions for denovo")

  class Params(define.Params):
  	model_id = define.string(required=True,
                             description='A unique name for this instance of the insert generator', category='General')
  	chromosome = define.integer(default=1, list=True)
  	p_del = define.real(default=0.01, min=0, max=1, category='Model params', description='probability of a deletion')
  	phet = define.real(default=0.05, min=0, max=1, category='Model params',
  	                   description='probability of having heterozygous mutation')
  	del_len_lo = define.integer(default=100, min=1, category='Model params',
  	                         description='Lower bound on deletion lengths')
  	del_len_high = define.integer(default=10000, min=1, category='Model params',
  	                         description='Upper bound on deletion lengths')
  	del_loc_rng_seed = define.integer(default=1, min=0, max=2**32 - 1, category='Model params: RNG',
  	      description='Seed for data location random number generator')
  	del_len_rng_seed = define.integer(default=2, min=0, max=2**32 - 1, category='Model params: RNG',
  	      description='Seed for deletion length random number generator')
  	het_rng_seed = define.integer(default=3, min=0, max=2**32 - 1, category='Model params: RNG',
  	      description='Seed for deletion length random number generator')
  	copy_rng_seed = define.integer(default=2, min=0, max=2**32 - 1, category='Model params: RNG')

  def write_to_json(self, fname):
  	with open(fname, 'w') as f:
  	  params = self.params.__json__()
  	  params.pop('model_id')
  	  json.dump({'delete': dict(**params)}, f)

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
  	"chromosome": [1],
  	"model_id": "delete_test",
    "phet": 0,
    "p_del": 0.01,
    "del_len_lo": 100,
    "del_len_high": 10000,
    "del_loc_rng_seed": 1,
    "del_len_rng_seed": 2,
    "het_rng_seed": 3,
    "copy_rng_seed": 2
  }
  inputs = {}
  wrp = Deletion(inputs, params)
  outputs = wrp.test()
  expected = {"delete": {
  		"chromosome": [1],
	    "phet": 0,
	    "p_del": 0.01,
	    "del_len_lo": 100,
	    "del_len_high": 10000,
	    "del_loc_rng_seed": 1,
	    "del_len_rng_seed": 2,
	    "het_rng_seed": 3,
	    "copy_rng_seed": 2
	    }
  }
  assert_equals.__self__.maxDiff = None
  with open(outputs.json_fragment) as fp:
    assert_equals(json.load(fp), expected)
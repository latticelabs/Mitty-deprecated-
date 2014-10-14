"""
This wrapper converts the parameters into a json file that is fed into the denovo wrap.

Schema example:

{
  "chromosome": [1],      # List of chromosomes to apply the variant to
  "p": 0.001,             # probability that the SNP will happen at any given base
  "phet": 0.5,            # probability that the variant will be heterozygous
  "base_loc_rng_seed": 1, # Seeds for the RNGs
  "base_sub_rng_seed": 2,
  "het_rng_seed": 3,
  "copy_rng_seed": 4
}
"""

import os
import json
from sbgsdk import define
from nose.tools import assert_equals

class SNP(define.Wrapper):
  class Inputs(define.Inputs):
    pass

  class Outputs(define.Outputs):
    json_fragment = define.output(name="SNP", description="SNP definitions for denovo")

  class Params(define.Params):
  	model_id = define.string(required=True,
                             description='A unique name for this instance of the insert generator', category='General')
  	chromosome = define.integer(default=1, list=True)
  	p = define.real(default=0.001, min=0, max=1, category='Model params', description='probability that the SNP will happen at any given base')
  	phet = define.real(default=0.05, min=0, max=1, category='Model params',
  	                   description='probability of having heterozygous mutation')
  	base_loc_rng_seed = define.integer(default=1, min=0, max=2**32 - 1, category='Model params: RNG',
  	      description='')
  	base_sub_rng_seed = define.integer(default=2, min=0, max=2**32 - 1, category='Model params: RNG',
  	      description='')
  	het_rng_seed = define.integer(default=3, min=0, max=2**32 - 1, category='Model params: RNG',
  	      description='Seed for deletion length random number generator')
  	copy_rng_seed = define.integer(default=2, min=0, max=2**32 - 1, category='Model params: RNG')

  def write_to_json(self, fname):
    with open(fname, 'w') as f:
      params = self.params.__json__()
      params.pop('model_id')
      json.dump({'snp': dict(**params)}, f)

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
    "chromosome": [1],      
    "p": 0.001,             
    "phet": 0.5,            
    "base_loc_rng_seed": 1, 
    "base_sub_rng_seed": 2,
    "het_rng_seed": 3,
    "copy_rng_seed": 4
  }
  inputs = {}
  wrp = SNP(inputs, params)
  outputs = wrp.test()
  expected = {"snp": {
      "chromosome": [1],    
      "p": 0.001,             
      "phet": 0.5,            
      "base_loc_rng_seed": 1, 
      "base_sub_rng_seed": 2,
      "het_rng_seed": 3,
      "copy_rng_seed": 4
	    }
  }
  assert_equals.__self__.maxDiff = None
  with open(outputs.json_fragment) as fp:
    assert_equals(json.load(fp), expected)
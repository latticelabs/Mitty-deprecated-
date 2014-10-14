"""
This wrapper converts the parameters into a json file that is fed into the denovo wrap.

Schema example:

{
  "chromosome": [1],
  "phet": 0.5,
  "p": 0.001,
  "ins_len_lo": 5,
  "ins_len_hi": 10
}
"""

import os
import json
from sbgsdk import define
from nose.tools import assert_equals

class Insert(define.Wrapper):
  class Inputs(define.Inputs):
    pass

  class Outputs(define.Outputs):
    json_fragment = define.output(name="Insertions", description="Insertion definitions for denovo")

  class Params(define.Params):
  	model_id = define.string(required=True,
                             description='A unique name for this instance of the insert generator', category='General')
  	chromosome = define.integer(default=1, list=True)
  	p = define.real(default=0.01, min=0, max=1, category='Model params', description='probability of an insertion')
  	phet = define.real(default=0.05, min=0, max=1, category='Model params',
  	                   description='probability of having heterozygous mutation')
  	ins_len_lo = define.integer(default=5, min=1, category='Model params',
  	                         description='Lower bound on insertion lengths')
  	ins_len_hi = define.integer(default=10, min=1, category='Model params',
  	                         description='Upper bound on insertion lengths')

  def write_to_json(self, fname):
  	with open(fname, 'w') as f:
  	  params = self.params.__json__()
  	  params.pop('model_id')
  	  json.dump({'insert': dict(**params)}, f)

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
    "chromosome": [1],
    "phet": 0.5,
    "p": 0.001,
    "ins_len_lo": 5,
    "ins_len_hi": 10,
    "model_id": "insert_test"
  }
  inputs = {}
  wrp = Insert(inputs, params)
  outputs = wrp.test()
  expected = {"insert": {
    "chromosome": [1],
    "phet": 0.5,
    "p": 0.001,
    "ins_len_lo": 5,
    "ins_len_hi": 10
	    }
  }
  assert_equals.__self__.maxDiff = None
  with open(outputs.json_fragment) as fp:
    assert_equals(json.load(fp), expected)
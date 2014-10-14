"""
This wrapper converts the parameters into a json file that is fed into the vcf2reads wrap.

Schema example:

{
  'paired': True,      #Are the reads paired or not.
  'read_len': 100,     #length of each read
  'template_len': 250, #length of template (only used for paired reads)
  'read_advance': 20   #how much to advance along the reference after generating a read. Determines "coverage"
}


"""

import os
import json
from sbgsdk import define
from nose.tools import assert_equals

class Simple_Seq(define.Wrapper):
  class Inputs(define.Inputs):
    pass

  class Outputs(define.Outputs):
    json_fragment = define.output(name="simple_sequential", description="Simple_seq definitions for vcf2_reads")

  class Params(define.Params):
    #model_id = define.string(required=True, description='A unique name', category='General')
    paired = define.boolean(name='paired reads', default='True')
    read_len = define.integer(name='read length', default=100)
    template_len = define.integer(name='length of template', description='only used for paired reads', default=250)
    read_advance = define.integer(name='read advance', description='how much to advance along the reference after generating a read. Determines coverage', default=20)


  def write_to_json(self, fname):
    with open(fname, 'w') as f:
      params = self.params.__json__()
      json.dump(dict(**params), f)

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    self.outputs.json_fragment = \
      os.path.join(output_dir, 'simple_sequential_plugin_params.json')
    self.write_to_json(self.outputs.json_fragment)
    self.outputs.json_fragment.meta = self.outputs.json_fragment.make_metadata(file_type='json')

def test():
  params = {
    "paired": True,
    "read_len": 100,
    "template_len": 250,
    "read_advance": 25
  }
  inputs = {}
  wrp = Simple_Seq(inputs, params)
  outputs = wrp.test()
  expected = {
  'paired': True,
  'read_len': 100,
  'template_len': 250,
  'read_advance': 25
    }

  assert_equals.__self__.maxDiff = None
  with open(outputs.json_fragment) as fp:
    assert_equals(json.load(fp), expected)

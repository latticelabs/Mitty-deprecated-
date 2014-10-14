"""
This wrapper converts the parameters into a json file that is fed into the vcf2reads wrap.

Schema example:

{
  'read_len': 100,     #length of each read
  'template_len': 250, #length of template
  'coverage': 1.0,     #How many x coverage do we want
  'max_p_error': 0.01,  #Maximum error rate at tip of read
  'k': 0.3
}


"""

import os
import json
from sbgsdk import define
from nose.tools import assert_equals

class Simple_Illumina(define.Wrapper):
  class Inputs(define.Inputs):
    pass

  class Outputs(define.Outputs):
    json_fragment = define.output(name="simple_illumina", description="Simple_illumina definitions for vcf2_reads")

  class Params(define.Params):
    read_len = define.integer(name='read length', default=100)
    coverage = define.real(name='coverage', default=1.0)
    template_len = define.integer(name='length of template', description='only used for paired reads', default=250)
    max_p_error = define.real(name='max_p_error', description='Maximum error rate at tip of read', default=0.01)
    k = define.real(name='k', default=0.3)

  def write_to_json(self, fname):
    with open(fname, 'w') as f:
      params = self.params.__json__()
      json.dump(dict(**params), f)

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    self.outputs.json_fragment = \
      os.path.join(output_dir, 'simple_illumina_plugin_params.json')
    self.write_to_json(self.outputs.json_fragment)
    self.outputs.json_fragment.meta = self.outputs.json_fragment.make_metadata(file_type='json')

def test():
  params = {
    "coverage": 1.0,
    "read_len": 100,
    "template_len": 250,
    "max_p_error": 0.01,
    "k": 0.3
  }
  inputs = {}
  wrp = Simple_Illumina(inputs, params)
  outputs = wrp.test()
  expected = {
    "coverage": 1.0,
    "read_len": 100,
    "template_len": 250,
    "max_p_error": 0.01,
    "k": 0.3
    }

  assert_equals.__self__.maxDiff = None
  with open(outputs.json_fragment) as fp:
    assert_equals(json.load(fp), expected)

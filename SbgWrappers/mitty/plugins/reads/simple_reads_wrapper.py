"""The snp_plugin does not run by itself: it is used to generate parts of the .json file that will be used by
reads_wrapper to create a full .json file that can be fed to reads.py

The .json fragment should look like

{
    "model": "simple_reads",
    "args": {
        "paired": false,
        "read_len": 100,
        "template_len": 250,
        "read_loc_rng_seed": 0,
        "error_rng_seed": 1,
        "base_chose_rng_seed": 2,
        "max_p_error": 0.8,
        "k": 0.1
    }
}

reads_wrapper will append this dictionary to the rest of the parameter file it is generating
"""
import os
import json
from sbgsdk import define
from nose.tools import assert_equals


class SimpleReads(define.Wrapper):
  # class Inputs(define.Inputs):
  #   pass

  class Outputs(define.Outputs):
    json_fragment = define.output(name="Simple Reads", description="Simple fixed length, random reads with exponential error profile")

  class Params(define.Params):
    paired = define.boolean(description='Are these paired end reads?', default=False, category='Template')
    read_len = define.integer(description='Length of reads', default=100, category='Template')
    template_len = define.integer(description='Length of template', default=250, category='Template')
    read_loc_rng_seed = define.integer(default=1, description='Seed for RNG governing where reads are taken from', required=True,
                                       min=0, max=2**32-1, category='RNG seeds')
    error_rng_seed = define.integer(default=2, description='Seed for RNG governing location of read errors', required=True,
                                    min=0, max=2**32-1, category='RNG seeds')
    base_chose_rng_seed = define.integer(default=3, description='Seed for RNG governing base substitution of read errors', required=True,
                                         min=0, max=2**32-1, category='RNG seeds')
    max_p_error = define.real(default=0.8,
                              description='Error level at tip end of read strand', required=True,
                              min=0, max=1.0, category='Read Errors')
    k = define.real(default=0.1, description='Exponential constant for read error envelope. The smaller this is, the faster the error rate vanishes to zero as we walk along the read from tip to base', required=True,
                              min=0, max=1.0, category='Read Errors')

  def write_to_json(self, fname):
    with open(fname, 'w') as f:
      params = self.params.__json__()
      json.dump({'read_model': 'simple_reads', 'model_params': params}, f)

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    self.outputs.json_fragment = \
      os.path.join(output_dir, 'read_plugin_params.json')  # Fixed name since we can't have multiple read models
    self.write_to_json(self.outputs.json_fragment)
    self.outputs.json_fragment.meta = self.outputs.json_fragment.make_metadata(file_type='json')


def test():
  params = {
    "paired": False,
    "read_len": 100,
    "template_len": 250,
    "read_loc_rng_seed": 0,
    "error_rng_seed": 1,
    "base_chose_rng_seed": 2,
    "max_p_error": 0.8,
    "k": 0.1
  }
  inputs = {}
  wrp = SimpleReads(inputs, params)
  outputs = wrp.test()
  expected = {
    "read_model": "simple_reads",
    "model_params": {
        "paired": False,
        "read_len": 100,
        "template_len": 250,
        "read_loc_rng_seed": 0,
        "error_rng_seed": 1,
        "base_chose_rng_seed": 2,
        "max_p_error": 0.8,
        "k": 0.1
    }
  }
  with open(outputs.json_fragment) as fp:
    assert_equals(json.load(fp), expected)

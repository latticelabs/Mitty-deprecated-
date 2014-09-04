"""Automatically find variant plugins and perform integration tests on them. For a plugin to avail of this test
the plugin needs a _example_params() function that returns a complete parameter file"""
import importlib
from nose.plugins.skip import SkipTest
import glob
import json
import mitty.denovo
from ... import *


def find_plugin_models():
  plugin_dir = glob.os.path.join(glob.os.path.dirname(__file__), glob.os.pardir, glob.os.pardir, glob.os.pardir, 'mitty', 'Plugins', 'variants')
  # The only awkward thing here is the truncation of '.py' ([:-3]) from the file name
  return [(mod, importlib.import_module('mitty.Plugins.variants.' + mod))
          for mod in [glob.os.path.basename(fn)[:-3] for fn in glob.glob(glob.os.path.join(plugin_dir, '*_plugin.py'))]]


def check_plugin_integration(model):
  if not hasattr(model[1], '_example_params'):
    #http://stackoverflow.com/questions/1120148/disabling-python-nosetests
    raise SkipTest('{:s} has no _example_params method. Can not test automatically'.format(model[0]))
  vcf_name = glob.os.path.join(data_dir, '{:s}_plugin_test.vcf.gz'.format(model[0]))
  param_file = glob.os.path.join(data_dir, '{:s}_plugin_test.json'.format(model[0]))
  params = model[1]._example_params()
  json.dump(params, open(param_file, 'w'), indent=2)
  mitty.denovo.main(wg_name, vcf_file_name=vcf_name, param_file_name=param_file, master_seed=0)
  assert os.path.exists(vcf_name)  # A very simple test to see if the plugin doesn't crash

  os.remove(vcf_name)  # Be neat
  os.remove(param_file)


#http://stackoverflow.com/questions/19071601/how-do-i-run-multiple-python-test-cases-in-a-loop
def integration_test_all_found_plugins():
  """Integration test on automatically found mutation plugin"""
  for model in find_plugin_models():
    check_plugin_integration.description = model[0] + ' integration test'
    yield check_plugin_integration, model


def check_plugin(model):
  if not hasattr(model[1], 'test'):
    #http://stackoverflow.com/questions/1120148/disabling-python-nosetests
    raise SkipTest('{:s} has no test method. Can not test automatically'.format(model[0]))
  model[1].test()


#http://stackoverflow.com/questions/19071601/how-do-i-run-multiple-python-test-cases-in-a-loop
def self_test_all_found_plugins():
  """Integration test on automatically found mutation plugin"""
  for model in find_plugin_models():
    check_plugin.description = model[0] + ' self test'
    yield check_plugin, model

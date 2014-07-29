"""Automatically find mutation plugins and perform integration tests on them. For a plugin to avail of this test
the plugin needs a _example_params() function that returns a complete parameter file"""
from nose.plugins.skip import SkipTest
import glob
import json
import mitty.mutate
from ... import *


def find_plugin_models():
  plugin_dir = glob.os.path.join(glob.os.path.dirname(__file__), glob.os.pardir, glob.os.pardir, glob.os.pardir, 'mitty', 'Plugins', 'Mutation')
  # The only awkard think here is the truncation of '_plugin.py' ([:-10]) from the file name
  return mitty.mutate.load_models([{'model': glob.os.path.basename(fn)[:-10]} for fn in glob.glob(glob.os.path.join(plugin_dir, '*_plugin.py'))])


def check_plugin(model):
  if not hasattr(model, '_example_params'):
    #http://stackoverflow.com/questions/1120148/disabling-python-nosetests
    raise SkipTest('{:s} has no _example_params method. Can not test automatically'.format(model.__name__))
  vcf_name = glob.os.path.join(data_dir , 'low_entropy_test.vcf.gz')
  param_file = glob.os.path.join(data_dir, 'low_entropy_test.json')
  params = model._example_params()
  json.dump(params, open(param_file, 'w'), indent=2)
  args = {
    '--wg': wg_name,
    '--vcf': vcf_name,
    '--paramfile': param_file,
    '--master_seed': 0,
    '-v': True
  }
  mitty.mutate.main(args)
  assert os.path.exists(vcf_name)  # A very simple test to see if the plugin doesn't crash


#http://stackoverflow.com/questions/19071601/how-do-i-run-multiple-python-test-cases-in-a-loop
def test_all_found_plugins():
  """Integration test on automatically found mutation plugin"""
  for model in find_plugin_models():
    check_plugin.description = 'Testing {:s} automatically'.format(model.__name__)
    yield check_plugin, model

"""Automatically find variant plugins and perform integration tests on them. For a plugin to avail of this test
the plugin needs a _example_params() function that returns a complete parameter file"""
from inspect import getmembers, isfunction
from nose.plugins.skip import SkipTest
import mitty.plugins.putil as putil
import mitty.lib.genome as genome
import glob
import json
import mitty.denovo
from ... import *


ref = genome.FastaGenome(example_fasta_genome)


def check_plugin_integration(model):
  if not hasattr(model[1], '_example_params'):
    #http://stackoverflow.com/questions/1120148/disabling-python-nosetests
    raise SkipTest('{:s} has no _example_params method. Can not test automatically'.format(model[0]))
  vcf_name = glob.os.path.join(data_dir, '{:s}_plugin_test.vcf.gz'.format(model[0]))
  params = {
    "denovo_variant_models": [
      {model[0]: model[1]._example_params}
    ]
  }
  mitty.denovo.main(ref, vcf_file_name=vcf_name, parameters=params, master_seed=1)
  assert os.path.exists(vcf_name)  # A very simple test to see if the plugin doesn't crash

  os.remove(vcf_name)  # Be neat


#http://stackoverflow.com/questions/19071601/how-do-i-run-multiple-python-test-cases-in-a-loop
def integration_test_all_found_plugins():
  """Integration test on automatically found mutation plugin"""
  for model in putil.load_all_variant_plugins():
    check_plugin_integration.description = model[0] + ' integration test'
    yield check_plugin_integration, model


def t_est_wrapper(func):
  func[1]()


def plugin_has_no_t_ests(model):
  raise SkipTest('No tests')


#http://stackoverflow.com/questions/19071601/how-do-i-run-multiple-python-test-cases-in-a-loop
def self_test_all_found_plugins():
  """Plugin self test"""
  for model in putil.load_all_variant_plugins():
    tests = [v for v in getmembers(model[1], isfunction) if v[0].startswith('test')]
    if len(tests) == 0:
      plugin_has_no_t_ests.description = model[0] + ' plugin self test(s)'
      yield plugin_has_no_t_ests, model
    else:
      for test in tests:
        t_est_wrapper.description = model[0] + ' plugin self test(s): ' + (test[1].func_doc or test[1].__name__)
        # We can't ensure that a dev will provide us with a function doc, so we use the name if can't find a doc string
        yield t_est_wrapper, test

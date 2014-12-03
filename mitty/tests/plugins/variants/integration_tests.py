"""Automatically find variant plugins and perform integration tests on them. For a plugin to avail of this test
the plugin needs a _example_params() function that returns a complete parameter file"""
from inspect import getmembers, isfunction
from nose.plugins.skip import SkipTest
from nose.tools import nottest
import mitty.lib.genome as genome
import mitty.denovo
from ... import *


ref = genome.FastaGenome(example_fasta_genome)


def check_plugin_integration(args):
  name, model = args
  if not hasattr(model, '_example_params'):
    #http://stackoverflow.com/questions/1120148/disabling-python-nosetests
    raise SkipTest('{:s} has no _example_params method. Can not test automatically'.format(name))
  params = {
    "denovo_variant_models": [
      {name: model._example_params}
    ]
  }
  g1 = mitty.denovo.main(ref, models=mitty.denovo.load_variant_model_list(params['denovo_variant_models']),
                         master_seed=1)
  assert type(g1) == dict  # A very simple test to see if the plugin doesn't crash


#http://stackoverflow.com/questions/19071601/how-do-i-run-multiple-python-test-cases-in-a-loop
def integration_test_all_found_plugins():
  """Integration test on automatically found mutation plugin"""
  for name, module in mitty.lib.discover_all_variant_plugins():
    check_plugin_integration.description = name + ' integration test'
    yield check_plugin_integration, (name, mitty.lib.load_variant_plugin(name))


@nottest
def test_wrapper(func):
  func[1]()


@nottest
def plugin_has_no_tests(_):
  raise SkipTest('No tests')


#http://stackoverflow.com/questions/19071601/how-do-i-run-multiple-python-test-cases-in-a-loop
def self_test_all_found_plugins():
  """Plugin self test"""
  for name, module in mitty.lib.discover_all_variant_plugins():
    model = mitty.lib.load_variant_plugin(name)
    tests = [v for v in getmembers(model, isfunction) if v[0].startswith('test')]
    if len(tests) == 0:
      plugin_has_no_tests.description = name + ' plugin self test(s)'
      yield plugin_has_no_tests, None
    else:
      for test in tests:
        test_wrapper.description = name + ' plugin self test(s): ' + (test[1].func_doc or test[1].__name__)
        # We can't ensure that a dev will provide us with a function doc, so we use the name if can't find a doc string
        yield test_wrapper, test

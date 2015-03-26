"""Automatically find variant plugins and perform integration tests on them. For a plugin to avail of this test
the plugin needs a _example_params() function that returns a complete parameter file"""
import tempfile
from inspect import getmembers, isfunction

from nose.plugins.skip import SkipTest
from nose.tools import nottest

import mitty.lib.db as mdb
import mitty.lib.io as mio
import mitty.genomes as genomes
from mitty.plugins.site_frequency import double_exp

from mitty.tests import *


ref = mio.Fasta(multi_dir=example_fasta_genome)


def check_plugin_integration(args):
  name, model = args
  if not hasattr(model, '_example_params'):
    #http://stackoverflow.com/questions/1120148/disabling-python-nosetests
    raise SkipTest('{:s} has no _example_params method. Can not test automatically'.format(name))
  params = [{name: model._example_params}]

  temp_fp, temp_name = tempfile.mkstemp(suffix='.sqlite3')
  os.close(temp_fp)
  pop_db_name = temp_name
  sfs_model = double_exp.Model()
  variant_models = genomes.load_variant_models(ref, params)
  chromosomes = [1]
  sample_size = 2
  master_seed = 2
  genomes.run_simulations(pop_db_name, ref, sfs_model, variant_models, chromosomes, sample_size, master_seed)

  # If we get here, the simulation ran. We just want a superficial test to round things out
  conn = mdb.connect(temp_name)
  ml = mdb.load_master_list(conn, 1)
  assert ml is not None


#http://stackoverflow.com/questions/19071601/how-do-i-run-multiple-python-test-cases-in-a-loop
def integration_test_all_found_plugins():
  """Integration test on automatically found mutation plugin"""
  for name, module in mitty.lib.discover_all_variant_plugins():
    check_plugin_integration.description = name + ' (variant plugin) integration test'
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
      plugin_has_no_tests.description = name + ' (variant plugin) self test(s)'
      yield plugin_has_no_tests, None
    else:
      for test in tests:
        test_wrapper.description = name + ' (variant plugin) self test(s): ' + (test[1].func_doc or test[1].__name__)
        # We can't ensure that a dev will provide us with a function doc, so we use the name if can't find a doc string
        yield test_wrapper, test

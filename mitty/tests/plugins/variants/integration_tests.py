"""Automatically find variant plugins and perform integration tests on them. For a plugin to avail of this test
the plugin needs a _example_params() function that returns a complete parameter file"""
import tempfile
import json
from inspect import getmembers, isfunction

from click.testing import CliRunner
from nose.plugins.skip import SkipTest
from nose.tools import nottest

import mitty.lib.variants as vr
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
  var_model_params = [{name: model._example_params}]

  _, param_file = tempfile.mkstemp(suffix='.json')
  _, db_file = tempfile.mkstemp(suffix='.hdf5')
  test_params = {
    "files": {
      "reference_dir": example_data_dir,
      "dbfile": db_file
    },
    "rng": {
      "master_seed": 1
    },
    "sample_size": 2,
    "chromosomes": [1, 2],
    "variant_models": var_model_params
  }
  json.dump(test_params, open(param_file, 'w'))

  runner = CliRunner()
  result = runner.invoke(genomes.cli, ['generate', param_file])
  assert result.exit_code == 0, result
  assert os.path.exists(db_file)

  pop = vr.Population(fname=db_file)
  ml = pop.get_master_list(chrom=1)

  os.remove(param_file)
  os.remove(db_file)

  assert len(ml) > 0


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


def sanity_check_all_found_plugins_test():
  """Sanity check on automatically found mutation plugin"""
  for name, module in mitty.lib.discover_all_variant_plugins():
    variant_sanity_check.description = name + ' (variant plugin) sanity check'
    yield variant_sanity_check, mitty.lib.load_variant_plugin(name).Model()


def variant_sanity_check(m):
  """Convenience function. Given an initialized model try and do a sanity check test with it."""
  ref_seq = ref[1]['seq']
  pos, stop, refs, alts, p = m.get_variants(ref_seq, seed=10)
  if len(pos) == 0:
    raise SkipTest('The defaults do not yield any variants to test')

  for p, s, r, a in zip(pos, stop, refs, alts):
    assert r[0] == ref_seq[p]
    if len(r) != len(a):
      assert a[0] == ref_seq[p]
    assert s == p + len(r)
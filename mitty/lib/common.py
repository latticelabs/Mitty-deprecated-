"""Common methods and constants imported when we import mitty.lib."""

import os
import random
import sys
import warnings
import pkg_resources
import inspect as py_inspect
from itertools import izip_longest  # Both for inspecting default values

SEED_MAX = (1 << 32) - 1  # Used for seeding rng
SFS_PLUGIN_ENTRY_POINT = 'mitty.plugins.sfs'
VARIANT_PLUGIN_ENTRY_POINT = 'mitty.plugins.variants'
POP_PLUGIN_ENTRY_POINT = 'mitty.plugins.population'
READS_PLUGIN_ENTRY_POINT = 'mitty.plugins.reads'
BENCHMARK_TOOL_WRAPPER_ENTRY_POINT = 'mitty.benchmarking.tools'


def rpath(base_dir, this_path):
  """Return this_path relative to base_dir, unless this_path is absolute"""
  if this_path is not None:
    return this_path if os.path.isabs(this_path) else os.path.expanduser(os.path.normpath(os.path.join(base_dir, this_path)))
  else:
    return None


def get_seeds(master_seed=1, size=1):
  """Use stock ``random`` to create `size` seeds from the master_seed."""
  r = random.Random(master_seed)
  return [r.randint(1, SEED_MAX) for _ in range(size)]


def discover_all_sfs_plugins():
  return sorted([(v.name, v.module_name) for v in pkg_resources.iter_entry_points(SFS_PLUGIN_ENTRY_POINT)],
                cmp=lambda x, y: cmp(x[0], y[0]))


def discover_all_variant_plugins():
  return sorted([(v.name, v.module_name) for v in pkg_resources.iter_entry_points(VARIANT_PLUGIN_ENTRY_POINT)],
                cmp=lambda x, y: cmp(x[0], y[0]))


def discover_all_reads_plugins():
  return sorted([(v.name, v.module_name) for v in pkg_resources.iter_entry_points(READS_PLUGIN_ENTRY_POINT)],
                cmp=lambda x, y: cmp(x[0], y[0]))


def _load_plugin(name, plugin_entry_point):
  v = [v1 for v1 in pkg_resources.iter_entry_points(plugin_entry_point, name)]
  if len(v) == 0:
    raise ImportError('No plugin called "{:s}" has been registered.'.format(name))
  if len(v) > 1:
    warnings.warn('More than one model with that name found. Loading first one only.')
  return v[0].load()


def load_sfs_plugin(name):
  return _load_plugin(name, SFS_PLUGIN_ENTRY_POINT)


def load_variant_plugin(name):
  return _load_plugin(name, VARIANT_PLUGIN_ENTRY_POINT)


def load_pop_model_plugin(name):
  return _load_plugin(name, POP_PLUGIN_ENTRY_POINT)


def load_reads_plugin(name):
  return _load_plugin(name, READS_PLUGIN_ENTRY_POINT)


def load_benchmark_tool_wrapper(name):
  for v in pkg_resources.iter_entry_points(BENCHMARK_TOOL_WRAPPER_ENTRY_POINT, name):
    return v.load()
  raise ImportError('No tool wrapper called "{:s}" has been registered.'.format(name))


import string
DNA_complement = string.maketrans('ATCGN', 'TAGCN')


def progress_bar(title, f, cols):
  """Draw a nifty progress bar.
  '\r' trick from http://stackoverflow.com/questions/15685063/print-a-progress-bar-processing-in-python

  :param title: leading text to print
  :param f:     fraction completed
  :param cols:  how many columns wide should the bar be
  """
  x = int(f * cols + 0.5)
  sys.stdout.write('\r' + title + '[' + '.' * x + ' ' * (cols - x) + ']')
  sys.stdout.flush()


def model_init_signature_string(obj):
  """Given an Object, show it's __init__ signature so we can get the parameter defaults
  :param obj: object
  """
  arg_spec = py_inspect.getargspec(obj)
  return 'Model defaults: (' + ',  '.join(reversed(['{:s}={:s}'.format(k, str(v)) for k, v in izip_longest(reversed(arg_spec.args), reversed(arg_spec.defaults), fillvalue='?') if k != 'self'])) + ')'

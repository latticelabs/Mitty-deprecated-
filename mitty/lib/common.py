"""Common methods and constants imported when we import mitty.lib."""

import os
import random
import warnings
import pkg_resources

SEED_MAX = (1 << 32) - 1  # Used for seeding rng
SFS_PLUGIN_ENTRY_POINT = 'mitty.plugins.sfs'
VARIANT_PLUGIN_ENTRY_POINT = 'mitty.plugins.variants'
READS_PLUGIN_ENTRY_POINT = 'mitty.plugins.reads'
BENCHMARK_TOOL_WRAPPER_ENTRY_POINT = 'mitty.benchmarking.tools'


def rpath(base_dir, this_path):
  """Return this_path relative to base_dir, unless this_path is absolute"""
  if this_path is not None:
    return this_path if os.path.isabs(this_path) else os.path.normpath(os.path.join(base_dir, this_path))
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


def load_reads_plugin(name):
  return _load_plugin(name, READS_PLUGIN_ENTRY_POINT)


def load_benchmark_tool_wrapper(name):
  for v in pkg_resources.iter_entry_points(BENCHMARK_TOOL_WRAPPER_ENTRY_POINT, name):
    return v.load()
  raise ImportError('No tool wrapper called "{:s}" has been registered.'.format(name))


import string
DNA_complement = string.maketrans('ATCGN', 'TAGCN')
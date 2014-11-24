import os
import random
import pkg_resources

SEED_MAX = (1 << 32) - 1  # Used for seeding rng
VARIANT_PLUGIN_ENTRY_POINT = 'mitty.plugins.variants'
READS_PLUGIN_ENTRY_POINT = 'mitty.plugins.reads'


def rpath(base_dir, this_path):
  """Return this_path relative to base_dir, unless this_path is absolute"""
  return this_path if os.path.isabs(this_path) else os.path.normpath(os.path.join(base_dir, this_path))


def get_seeds(master_seed=1, size=1):
  """Use stock ``random`` to create `size` seeds from the master_seed."""
  r = random.Random(master_seed)
  return [r.randint(1, SEED_MAX) for _ in range(size)]


def discover_all_variant_plugins():
  return sorted([(v.name, v.module_name) for v in pkg_resources.iter_entry_points(VARIANT_PLUGIN_ENTRY_POINT)],
                cmp=lambda x, y: cmp(x[0], y[0]))


def discover_all_reads_plugins():
  return sorted([(v.name, v.module_name) for v in pkg_resources.iter_entry_points(READS_PLUGIN_ENTRY_POINT)],
                cmp=lambda x, y: cmp(x[0], y[0]))


def _load_plugin(name, plugin_entry_point):
  for v in pkg_resources.iter_entry_points(plugin_entry_point):
    if v.name == name:
      return v.load()
  raise ImportError('No plugin called "{:s}" has been registered.'.format(name))


def load_variant_plugin(name):
  return _load_plugin(name, VARIANT_PLUGIN_ENTRY_POINT)


def load_reads_plugin(name):
  return _load_plugin(name, READS_PLUGIN_ENTRY_POINT)
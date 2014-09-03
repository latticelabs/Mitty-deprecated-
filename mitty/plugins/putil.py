"""Useful utility functions for dealing with plugins."""
import pkgutil
import importlib
import os
variant_plugin_dir = os.path.join(os.path.dirname(__file__), 'variants')
read_plugin_dir = os.path.join(os.path.dirname(__file__), 'reads')


def load_variant_plugin(plugin_name):
  return importlib.import_module('mitty.plugins.variants.' + plugin_name + '_plugin')


def list_all_variant_plugins():
  return [package_name[:-7]
          for _, package_name, _ in pkgutil.iter_modules([variant_plugin_dir]) if package_name.endswith('_plugin')]


def load_all_variant_plugins():
  return [(name, load_variant_plugin(name)) for name in list_all_variant_plugins()]


def load_read_plugin(plugin_name):
  return importlib.import_module('mitty.plugins.reads.' + plugin_name + '_plugin')


def list_all_read_plugins():
  return [package_name[:-7]
          for _, package_name, _ in pkgutil.iter_modules([read_plugin_dir]) if package_name.endswith('_plugin')]


def load_all_read_plugins():
  return [(name, load_read_plugin(name)) for name in list_all_read_plugins()]

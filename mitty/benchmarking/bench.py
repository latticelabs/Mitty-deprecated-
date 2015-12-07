"""This module provides high level functions to control benchmarking tasks on the platform

An example benchmarking setup and run script is (also found in examples/benchmarking/ex1.py)

metadata = {
  "bench_run": "R1",
  "bench_name": "B1",
  "bench_inputs": {  -> ordered dictionary of inputs that vary to create benchmark runs
    <input_name1>: <file_info>,
    <input_name2>: <file_info>,
    <input_name3>: <file_info>
  },
  "tool": "gral-0.1.0"
}

file_info = {
  "tag": <short name we use>
  .... any metadata we want ....
  .... don't use ordered dicts here ....
}


"""
import json
import logging
import os
from collections import OrderedDict
from copy import deepcopy
import hashlib

logger = logging.getLogger(__name__)


def meta2tuplelist(meta):
  """We need to make sure that we save our Ordered Dicts properly. We convert them into list tuples just before
  we need to save them as json

  :param meta:
  :return:
  """
  return [(k, v.items()) if type(v) == OrderedDict else (k, v) for k, v in meta.items()]


def tuplelist2meta(json_meta):
  """Given a dict from meta2dict saved as a json then loaded, reconstruct the ordered dicts properly

  :param json_meta: (This will be a list of tuples, as created by meta2tuplelist
  :return:
  """
  return OrderedDict([(k, OrderedDict(v)) if k == 'inputs' else (k, v) for k, v in json_meta])


def create_filename_prefix_from_metadata(meta, use_hash=True):
  """Given a metadata dictionary, create a filename from it
  :param meta: benchmark metadata dictionary.
  :param use_hash: If True, instead of using a long filename, we will use the md5 hash of the filename instead

  <bench_run>.<bench_name>.<input_name1>-<file_tag>. (repeated as needed) .<tool>.<ext>
  """
  def flatten_input_meta(k, v):
    return [_k + '-' + _v['tag'] for _k, _v in v.items()] if k == 'bench_inputs' else [v]

  human_readable_string = '.'.join([x for k1, v1 in meta.items() for x in flatten_input_meta(k1, v1)])
  return hashlib.md5(human_readable_string).hexdigest() if use_hash else human_readable_string


def create_bench_spec(name, description,
                      file_list,
                      bench_combinations,
                      benchmark_tools,
                      tool_output_suffix):
  """Given some descriptions, place them into a package we call a benchmark spec

  :param name:
  :param description:
  :param file_list:
  :param bench_combinations:
  :param bench_mark_tools:
  :param tool_output_suffix:
  :return: a dict

  Principally, we do this to ensure that the "combinations" dict is an ordered dicts.
  It is important to keep the ordering consistent in order for the metadata -> filename mapping
  to work consistently."""
  return {
    "bench_name": name,
    "description": description,
    "file_list": {k: dict(v.items() + [('tag', k)]) for k, v in file_list.items()},
    "combinations": bench_combinations,
    "benchmark_tools": benchmark_tools,
    "tool_output_suffix": tool_output_suffix,
  }


def create_bench_run(name, description, bench_spec,
                     tool_descriptions, previous_bench_run=[],
                     use_hash_for_filenames=True):
  """Add tool descriptions to a bench spec and create all the run combinations required. The information
  should be sufficient to allow a task manager to create tasks for these.

  :param name:
  :param description:
  :param bench_spec:
  :param tool_descriptions:
  :param previous_bench_run:
  :return:
  """
  bench_run_spec = deepcopy(bench_spec)
  bench_run_spec['bench_run_name'] = name
  bench_run_spec['description'] = description
  bench_run_spec['tool_descriptions'] = OrderedDict([(k, dict(v.items() + [('tag', k)]))
                                                     for k, v in tool_descriptions])
  bench_run_spec['previous_run'] = previous_bench_run  # TODO: Implement this

  bench_run_spec['tool_and_analysis_task_list'] = \
    {k: tool_and_analysis_task_list(bench_run_spec, td, use_hash_for_filenames)
     for k, td in bench_run_spec['tool_descriptions'].items()}

  bench_run_spec['meta_analysis_task'] = compute_meta_analysis_task(bench_run_spec, use_hash_for_filenames)

  return bench_run_spec


def tool_and_analysis_task_list(bench_run_spec, tool_description, use_hash):
  """Given a bench_spec and a tool description, construct the tool and tool_analysis task descriptions"""
  # First go through the bench combinations and decide which entries are needed, then create tasks based on
  # only those
  return [compute_tool_and_analysis_task(bench_run_spec, td, tool_description, use_hash)
          for td in compute_all_combinations([(k, v) for k, v in bench_run_spec['combinations']
                                              if k in tool_description['input_mapping']])]


def compute_all_combinations(spec_list, this_run=OrderedDict()):
  """A recursive function to compute all combinations of lists found in inputs and tools

  :param spec_list: A list of tuples
  :param this_run: a dictionary of inputs and tool id
  :return: A list
  """
  def next_combination(_key, _val, _run):
    _new_run = deepcopy(_run)
    _new_run[_key] = _val
    return _new_run

  # Do the base case first - easier to reason about
  if len(spec_list) == 0:
    return [this_run]

  l = []
  for val in spec_list[0][1]:
    l += compute_all_combinations(spec_list[1:], next_combination(spec_list[0][0], val, this_run))
  return l


def compute_tool_and_analysis_task(bench_run_spec, task_dict, tool_desc, use_hash):
  tool_task = compute_tool_task(bench_run_spec, task_dict, tool_desc, use_hash)
  anal_task = compute_analysis_task(bench_run_spec, tool_desc, tool_task, use_hash)
  return {
    'tool_task': tool_task,
    'anal_task': anal_task
  }


def compute_tool_task(bench_run_spec, task_dict, tool_desc, use_hash):
  _files = {k: bench_run_spec['file_list'][v] for k, v in task_dict.items()}  # Inputs described in task_dict
  _files.update({k: bench_run_spec['file_list'][k] for k, _ in tool_desc['input_mapping'].items()
                 if k not in _files})  # Fixed inputs not part of the bench combinations
  input_files = {tool_desc['input_mapping'][k]: v for k, v in _files.items()}
  metadata = OrderedDict([
    ("bench_run", bench_run_spec['bench_run_name']),
    ("bench_name", bench_run_spec['bench_name']),
    ("bench_inputs", OrderedDict([(k, bench_run_spec['file_list'][v]) for k, v in task_dict.items()])),
    ("tool", tool_desc['tag'])
  ])
  output_files = {v: create_filename_prefix_from_metadata(metadata, use_hash) + '.' +
                     bench_run_spec['tool_output_suffix'][k]
                  for k, v in tool_desc['output_mapping'].items()}
  return {
    'input_files': input_files,
    'metadata': metadata,
    'output_files': output_files
  }


def compute_analysis_task(bench_run_spec, tool_description, tool_task, use_hash):
  """Inputs are picked from tool task outputs and input files"""
  fl = bench_run_spec['file_list']
  tom = tool_description['output_mapping']  # tool_output_mapping
  anal_inputs = bench_run_spec['benchmark_tools']['tool_analysis']['inputs']
  anal_outputs = bench_run_spec['benchmark_tools']['tool_analysis']['outputs']

  input_files = {k: tool_task['output_files'].get(tom.get(k, 'not a tool output'), fl.get(k, None))
                 for k in anal_inputs}
  metadata = deepcopy(tool_task['metadata'])
  output_files = {k: create_filename_prefix_from_metadata(metadata, use_hash) + '.' + v
                  for k, v in anal_outputs.items()}
  return {
    'input_files': input_files,
    'metadata': metadata,
    'output_files': output_files
  }


def compute_meta_analysis_task(bench_run_spec, use_hash):
  fl = bench_run_spec['file_list']
  tal = bench_run_spec['tool_and_analysis_task_list']
  mal_inputs = bench_run_spec['benchmark_tools']['meta_analysis']['inputs']
  mal_outputs = bench_run_spec['benchmark_tools']['meta_analysis']['outputs']

  input_files = {k: tv['tool_task']['output_files'].get(k, tv['anal_task']['output_files'].get(k, fl.get(k, None)))
                 for k in mal_inputs for tool, tool_tasks in tal.items() for tv in tool_tasks}
  metadata = deepcopy(next(tal.itervalues())[0]['tool_task']['metadata'])
  # This will fail if we have no tasks ...
  metadata['bench_inputs'] = {'many': {'tag': 'inputs'}}
  metadata['tool'] = 'many-tools'
  output_files = {k: create_filename_prefix_from_metadata(metadata, use_hash) + '.' + v
                  for k, v in mal_outputs.items()}
  return {
    'input_files': input_files,
    'metadata': metadata,
    'output_files': output_files
  }



# Future expansion - exclude list
def create_run_list(bench_run):
  """Compute all combinations of input files and tools to create a list of runs"""
  spec_list = convert_bench_spec_to_spec_list(bench_run['bench_spec']['input_files'])
  spec_list += [('app_to_bench', [k for k, _ in bench_run['tools'].iteritems()])]
  return compute_all_combinations(spec_list)




def create_tool_task_from_run(bench_run, this_run, use_hash=True):
  """Create explicit metadata and give extra details for the tool run

  :param bench_run: Benchmark details
  :param this_run: dictionary of input files and app.
                   'app_to_bench' is the reserved key used to indicate the app rather than the tools
  :return:
  """
  # The tool info
  tr = {
    'inputs': {},
    'metadata': {
      'benchmark_instance': bench_run['name'],
      'benchmark_name': bench_run['bench_spec']['name'],
      'inputs': OrderedDict([(k, {v: bench_run['bench_spec']['input_files'][k][v]}) for k, v in this_run.iteritems() if k != 'app_to_bench']),
      'tool': this_run['app_to_bench']
    }
  }

  out_suff = bench_run['bench_spec']['tool-output-suffix']
  # Figure out the output file names. Look through the output files in the tool description and find the suffix
  tr['outputs'] = {v: create_filename_prefix_from_metadata(tr['metadata'], use_hash=use_hash) + '.' + out_suff[k] for k, v in bench_run['tools'][this_run['app_to_bench']]['outputs'].iteritems()}

  for k, v in this_run.iteritems():
    if k is not 'app_to_bench':
      tr['inputs'][k] = bench_run['bench_spec']['input_files'][k][v]
      # This can't error unless we have made a mistake in chaining the commands together
    else:
      tr[k] = bench_run['tools'][v]

  # The analysis run
  ar = deepcopy(bench_run['bench_spec']['analysis'])
  for k in ['inputs', 'outputs', 'metadata']:
    if k in ar:
      RuntimeError('The analysis dictionary has a key "{}" which is reserved'.format(k))
  ar['inputs'] = {k: tr['outputs'][bench_run['tools'][this_run['app_to_bench']]['outputs'][k]]
                  for k in ar['tool-input-names']}
  ar['inputs'].update({k: bench_run['bench_spec']['other_files'][k]['file_name'] for k in ar['other-input-names']})

  out_suff = ar['output-suffix']
  ar['metadata'] = deepcopy(tr['metadata'])
  ar['metadata']['analysis'] = ar.get('name', '')
  ar['outputs'] = {k: create_filename_prefix_from_metadata(ar['metadata'], use_hash=use_hash) + '.' + v
                   for k, v in out_suff.iteritems()}

  return tr, ar


def create_meta_analysis_task(bench_run, platform_run_list, use_hash=True):
  """The meta analysis task list of inputs and the outputs

  :param bench_run:
  :param platform_run_list:
  :param use_hash:
  :return:
  """
  def find_file(_task, _inp_type):
    """Look through tool task and analysis task outputs to find a match for the file we need"""
    # This requires translation - it's the tool task and tools can have non-standard output mappings
    if _inp_type in _task[0]['app_to_bench']['outputs']:
      __inp_type = _task[0]['app_to_bench']['outputs'][_inp_type]
      return {'file-name': _task[0]['outputs'][__inp_type], 'metadata': _task[0]['metadata']}

    if _inp_type in _task[1]['outputs']:
      return {'file-name': _task[1]['outputs'][_inp_type], 'metadata': _task[1]['metadata']}

    return {'file-name': '', 'metadata': {}}

  mat = {
    'app_name': bench_run['bench_spec']['meta-analysis']['app_name'],
    'inputs': {inp_type: [find_file(task, inp_type) for task in platform_run_list]
               for inp_type in bench_run['bench_spec']['meta-analysis']['inputs']},
    'metadata': {
      'benchmark_instance': bench_run['name'],
      'benchmark_name': bench_run['bench_spec']['name'],
      'inputs': {'meta': {'many': None}},
      'tool': 'many',
      'analysis': bench_run['bench_spec']['analysis'].get('name', ''),
      'meta-analysis': bench_run['bench_spec']['meta-analysis'].get('name', '')
    }
  }
  mat['outputs'] = {k: create_filename_prefix_from_metadata(mat['metadata'], use_hash=use_hash) + '.' + v
                    for k, v in bench_run['bench_spec']['meta-analysis']['output-suffix'].iteritems()}
  return mat


def create_benchmarking_tasks(bench_run, run_list, use_hash=True):
  """Create a list of tool runs, analysis runs and meta-analysis runs with app names and input/output names as well
  as appropriate meta-data copied over. We could have bundled this with compute_all_combinations but it would have
  made the function more messy.

  :param bench_run:
  :param run_list: from create_run_list
  :param use_hash: Filenames use hashes instead of human readable names
  :return:
  """
  plat_tool_and_anal_task_list = [create_tool_task_from_run(bench_run, this_run, use_hash=use_hash) for this_run in run_list]
  plat_meta_anal_task = create_meta_analysis_task(bench_run, plat_tool_and_anal_task_list, use_hash=use_hash)

  return plat_tool_and_anal_task_list, plat_meta_anal_task

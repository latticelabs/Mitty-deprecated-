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

# These are for the BaseExecutor
import random
import time
from multiprocessing import Process


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

  # bench_run_spec['tool_and_analysis_task_list'] = OrderedDict(
  #   [(k, tool_and_analysis_task_list(bench_run_spec, td, use_hash_for_filenames))
  #    for k, td in bench_run_spec['tool_descriptions'].items()])

  bench_run_spec['tool_and_analysis_task_list'] = \
    [task for k, td in bench_run_spec['tool_descriptions'].items()
     for task in tool_and_analysis_task_list(bench_run_spec, td, use_hash_for_filenames)]

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
    'app': tool_desc,
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
    'app': bench_run_spec['benchmark_tools']['tool_analysis'],
    'input_files': input_files,
    'metadata': metadata,
    'output_files': output_files
  }


def compute_meta_analysis_task(bench_run_spec, use_hash):
  fl = bench_run_spec['file_list']
  tal = bench_run_spec['tool_and_analysis_task_list']
  mal_inputs = bench_run_spec['benchmark_tools']['meta_analysis']['inputs']
  mal_outputs = bench_run_spec['benchmark_tools']['meta_analysis']['outputs']

  # input_files = {k: tv['tool_task']['output_files'].get(k, tv['anal_task']['output_files'].get(k, fl.get(k, None)))
  #                for k in mal_inputs for tool, tool_tasks in tal.items() for tv in tool_tasks}
  input_files = {k: tv['tool_task']['output_files'].get(k, tv['anal_task']['output_files'].get(k, fl.get(k, None)))
                 for k in mal_inputs for tv in tal}

  #metadata = deepcopy(next(tal.itervalues())[0]['tool_task']['metadata'])
  metadata = deepcopy(tal[0]['tool_task']['metadata'])

  # This will fail if we have no tasks ...
  metadata['bench_inputs'] = {'many': {'tag': 'inputs'}}
  metadata['tool'] = 'many-tools'
  output_files = {k: create_filename_prefix_from_metadata(metadata, use_hash) + '.' + v
                  for k, v in mal_outputs.items()}
  return {
    'app': bench_run_spec['benchmark_tools']['meta_analysis'],
    'input_files': input_files,
    'metadata': metadata,
    'output_files': output_files
  }


# Benchmark controller -------------------------------------------------------------------------------------

def create_switch_board():
  """The dictionary that implements the state machine for the benchmark controller"""
  return {
    'global_state': {
      'waiting': start_tool_tasks,
      'tool_tasks': advance_tool_tasks,
      'meta_analysis_task': advance_meta_analysis_task,
      'done': nop,
      'error': nop
    }
  }


class BenchRunState(dict):
  """The only reason we do this is so that we can override the __repr__ method to print the state nicely"""
  def __repr__(self):
    return pprint_state(self)


def pprint_state(bench_run_state):
  """Pretty print the state.

  :param bench_run_state:
  :return:
  """
  lines = [
    "Bench run name: {bench_run_spec[bench_run_name]:s}",
    "Bench spec name: {bench_run_spec[bench_name]:s}",
    "Status: {global_state:s}",
    "-----------------------------------------------",
    "Tool runs: (Future versions will have short names for each run)",
    "----------------"
  ]
  lines += ["\tJob {job_id:s}: {state:s}".format(**task) for task in bench_run_state['task_states']['tool']]

  lines += [
    "Analysis runs:",
    "----------------"
  ]
  lines += ["\tJob {job_id:s}: {state:s}".format(**task) for task in bench_run_state['task_states']['analysis']]

  lines += ["Meta-analysis run: {task_states[meta-analysis][job_id]:s}: {task_states[meta-analysis][state]:s}"]

  return '\n'.join(lines).format(**bench_run_state)


def prepare_benchmark_run_state(bench_run_spec=None, old_bench_run_state=None):
  """Given a bench run spec, setup a state file that can be used to keep track of progress
  through the analysis

  :param bench_run_spec: The benchmark run spec file. If omitted, old_bench_run_state can not be None
  :param old_bench_run_state: If this is passed it uses this (as current state) instead of creating a new
                              state file (which starts from the default initial state)
                              Use this if loading a state from file
  :return: The benchmark run state
  """
  if old_bench_run_state is None:  # Create a fresh new state
    tal = bench_run_spec['tool_and_analysis_task_list']
    bench_run_state = {
      'bench_run_spec': bench_run_spec,
      'global_state': 'waiting',
      'task_states': {k: {'job_id': None, 'state': 'waiting'} if k == 'meta-analysis'
                          else [{'job_id': None, 'state': 'waiting'} for _ in tal]
                      for k in ['tool', 'analysis', 'meta-analysis']}
    }
  else:
    bench_run_state = deepcopy(old_bench_run_state)

  return BenchRunState(bench_run_state)


def bnext(bench_run_state, executor):
  """Given a state file, advance us to the next state. This is the only function we have to call

  :param bench_run_state: the state dictionary
  :param executor: the executor class that is called to run the tasks
  :return: the new state
  """
  switch_board = create_switch_board()
  return switch_board['global_state'][bench_run_state['global_state']](bench_run_state, executor)


def start_tool_tasks(bench_run_state, runner):
  bs = bench_run_state['bench_run_spec']
  tal = bs['tool_and_analysis_task_list']
  new_bench_run_state = deepcopy(bench_run_state)
  new_bench_run_state['task_states'].update({
    'tool': [{'job_id': runner.start_job(t['tool_task']['app'],
                                         t['tool_task']['input_files'],
                                         t['tool_task']['output_files']),
              'state': 'running'} for t in tal]})
  new_bench_run_state['global_state'] = 'tool_tasks'
  return new_bench_run_state


def advance_tool_tasks(bench_run_state, runner):
  """Given that some tool tasks are running, check on their status and if they are done, start the corresponding
  analysis tasks

  :param bench_run_state:
  :param runner:
  :return:
  """
  new_bench_run_state = deepcopy(bench_run_state)
  bs = new_bench_run_state['bench_run_spec']
  tt_list = new_bench_run_state['task_states']['tool']
  at_list = new_bench_run_state['task_states']['analysis']
  tal = bs['tool_and_analysis_task_list']

  #TODO: what to do with failed tasks?
  for tt, at, td in zip(tt_list, at_list, tal):
    tt['state'] = runner.job_status(tt['job_id'])
    if at['job_id'] is not None:
      at['state'] = runner.job_status(at['job_id'])
    if tt['state'] == 'finished' and at['state'] == 'waiting':
      at['job_id'] = runner.start_job(td['tool_task']['app'],
                                      td['tool_task']['input_files'],
                                      td['tool_task']['output_files'])
      at['state'] = 'running'

  # If all the analysis tasks are done, it's time to move on to meta-analysis
  if len(filter(lambda x: x['state'] in ['waiting', 'running'], at_list)) == 0:
    new_bench_run_state['global_state'] = 'meta_analysis_task'

  return new_bench_run_state


def advance_meta_analysis_task(bench_run_state, runner):
  """Given that all tool and analysis tasks have finished, start the meta-analysis task

  :param bench_run_state:
  :param runner:
  :return:
  """
  new_bench_run_state = deepcopy(bench_run_state)
  bs = new_bench_run_state['bench_run_spec']
  mat_stat = new_bench_run_state['task_states']['meta-analysis']
  mat = bs['meta_analysis_task']
  if mat_stat['state'] == 'waiting':
    mat_stat['job_id'] = runner.start_job(mat['app'],
                                          mat['input_files'],
                                          mat['output_files'])
    mat_stat['state'] = 'running'
  else:
    mat_stat['state'] = runner.job_status(mat_stat['job_id'])

  if mat_stat['state'] == 'finished':
    new_bench_run_state['global_state'] = 'done'

  return new_bench_run_state


def nop(bench_run_state, runner):
  """Do nothing, don't change state

  :param bench_run_state:
  :param runner:
  :return:
  """
  return bench_run_state


class BaseExecutor:
  """This is a dummy class which can be inherited to implement the function calls and state required by an
    actual runner. The default implementation simulates generating dummy files and can be used to test the
    framework."""

  def __init__(self, sleep_min=5, sleep_max=10):
    """Our fake tasks will take between sleep_min and sleep_max sec to run"""
    self.sleep_min, self.sleep_max = sleep_min, sleep_max
    self.job_id = {}  # Dict of processes indexed by PID

  def start_job(self, app, input_files, output_files):
    """This is where the

    :param app:
    :param input_files:
    :param output_files:
    :return:
    """
    p = Process(target=time.sleep, args=(random.randint(self.sleep_min, self.sleep_max),))
    p.start()
    pid = str(p.pid)
    self.job_id[pid] = p
    return pid

  def job_status(self, job_id):
    """

    :param job_id:
    :return: 'running', 'finished' or 'error'
    """
    if self.job_id[job_id].is_alive():
      return 'running'
    else:
      return 'finished'

if __name__ == '__main__':
  # Run a dummy task
  pass
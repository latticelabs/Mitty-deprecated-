"""This uses the 'Benchmarking-demo' project which has some dummy apps to help test the API interface."""
import mitty.benchmarking.bench as bench

file_list = {
  "i0": {
     "file_name": "input0.txt",
     "description": "Common file used by all tools"
  },
  "i1": {
   "file_name": "input1.txt",
   "meta1": "m11",
   "meta2": "m21",
   "description": "Input file with some metadata filled out"
  },
  "i2": {
   "file_name": "input2.txt",
   "meta1": "m12",
   "meta2": "m22",
   "description": "Input file with some metadata filled out"
  },
  "i3": {
   "file_name": "input3.txt",
   "metaA": "mA3",
   "metaB": "mA3",
   "description": "Input file with different metadata filled out"
  },
  "i4": {
   "file_name": "input4.txt",
   "metaA": "mA4",
   "metaB": "mA4",
   "description": "Input file with different metadata filled out"
  },
  "i5": {
   "file_name": "input5.txt",
   "description": "Input file only used by one tool"
  },
  "i6": {
   "file_name": "input6.txt",
   "description": "Input file only used by analysis"
  },
  "i7": {
   "file_name": "input7.txt",
   "description": "Input file only used by meta-analysis"
  }
}

"""
[ (key, [tag1, tag2, ....] .... )

Each tuple is an input dimension
key           - a string that maps onto a tool input
tag1,tag2,..  - list of file tags. Each tag MUST be found in file_list. For each dimension we loop over the
                files in this list
"""


bench_combinations = [
  ("iset1", ["i1", "i2"]),
  ("iset2", ["i3", "i4"])
]

benchmark_tools = {
  "tool_analysis": {
    # Description of the tool used to analyze each run. Leave out to indicate no per-tool analysis is done
    # (e.g. for comparative benchmarks where there is no truth set, we go directly to meta-analysis).
    "app_name": "analysis-tool",
    "tag": "anal-tool",
    "inputs": ["tool-out1", "tool-out2", "i6"],
    # "inputs" will be searched for in outputs of tools and in file_list
    # The name (e.g "bam") must be present in the "output_mapping" field of the tool description
    "outputs": {  # In general, we need to know the suffix for each file we produce.
      "anal-out1": "anal.out1.txt",
      "anal-out2": "anal.out2.txt"
    }
  },
  "meta_analysis": {
    # Description of the tool used to analyze all runs combined.
    # Leave out to skip this step
    "app_name": "meta-analysis-tool",
    "tag": "meta-anal-tool",
    "inputs": ["anal-out1", "anal-out2", "tool-out1", "i7"],
    # "inputs" will be searched for in outputs of tools, outputs of tool_analysis and in file_list
    # The name (e.g "tool-time") must be present in the "output_mapping" field of the tool description
    "outputs": {  # In general, we need to know the suffix for each file we produce.
      "meta-anal-out": "meta-out.txt"
    }
  }
}


# We rename all relevant tool outputs according to our naming convention.
# This mapping indicates what suffix we should use for relevant tool files
# It would be tedious to place this in the output description of each tool
# These are the outputs required and used by the benchmark tools
tool_output_suffix = {
  "tool-out1": "tool-out1.txt",
  "tool-out2": "tool-out2.txt"
}

bench_spec = bench.create_bench_spec(name='B1', description='Example SDK2 benchmark',
                                     file_list=file_list, bench_combinations=bench_combinations,
                                     benchmark_tools=benchmark_tools,
                                     tool_output_suffix=tool_output_suffix)


#TODO: How to handle additional files, like indexes, naturally. Ask Boysha
tool_descriptions = [
  ("tool1", {
    "app_name": "tool1",
    "tag": "t1",
    "input_mapping": {
      "i0": "tool1-input1",
      "iset1": "tool1-input2",
      "iset2": "tool1-input3"
    },
    "output_mapping": {
      "tool-out1": "tool1-output1",
      "tool-out2": "tool1-output2",
    },
    "parameters": {'-p1': 'p1'},  # Leave blank to leave defaults
    "description": "Tool 1"
  }),
  ("tool2", {
    "app_name": "tool2",
    "tag": "t2",
    "input_mapping": {
      "i0": "tool2-input1",
      "iset2": "tool2-input2"
    },
    "output_mapping": {
      "tool-out1": "tool2-output1",
      "tool-out2": "tool2-output2",
    },
    "parameters": {'-p1': 'p1'},  # Leave blank to leave defaults
    "description": "Tool 2"
  })
]

import json
import os
import tempfile
import logging
import mitty.benchmarking.sbgapi2 as api


d = tempfile.mkdtemp(prefix='temp_ex2', dir='./')  # Don't check this in BTW
os.chdir(d)
logging.basicConfig(filename='bench_test.log', level=logging.DEBUG)

bench_run = bench.create_bench_run(name='R1', description='Run1 with bench spec B1',
                                   bench_spec=bench_spec, tool_descriptions=tool_descriptions,
                                   use_hash_for_filenames=False)

json.dump(bench_run, open('bench_run.json', 'w'), indent=2)

exe = api.SBGSDK2Executor(bench_run, 'Bench Test')
#exe = bench.BaseExecutor(sleep_min=0.5, sleep_max=1)
bench.execute_benchmark(bench_run, exe, file_name='bench_state.json', poll_interval=0.1)

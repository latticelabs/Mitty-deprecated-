"""This sets up a benchmark and a benchmark run using some completely fake data files and the default
dummy runner. It generates a bunch of one line text files, so it's nice to run it in it's own directory
which removes the possibility it overwrites an actual data file!"""
import time

import mitty.benchmarking.bench as bench

# Let's set up our benchmark

# Describe the files we will have in the benchmark set. The format is a list of tuples which later get converted to
# an ordered dict. The keys used should match the names we use for the benchmarking

"""
A dict of all input files relevant to this benchmark. Each entry is of the form

tag: {
  "file_name": ... ,
  "description": ...,
  m1: ...,
  m2: ...,
  ...
}

tag (which serves as the key) is a short way to refer to the files and is referred to in the "inputs" section of
tool descriptions if the input should get this file directly.

description is a string meant for human consumption

m1, m2 are string keys whose values indicate metadata we need for the benchmarking
"""
file_list = {
  "ref": {
     "file_name": "ref.fa.gz",
     "description": "Dummy reference file"
  },
  "s0-r1": {
   "file_name": "reads.fq",
   "sample": "g0_s0",
   "read_type": "read-perfect",
   "description": "Perfect reads from sample g0_s0 across whole genome"
  },
  "s0-r2": {
   "file_name": "reads_c.fq",
   "sample": "v0",
   "read_type": "read-corrupt",
   "description": "Corrupt reads from sample g0_s0 across whole genome"
  },
  "s0-r3": {
   "file_name": "reads1.fq",
   "sample": "g0_s1",
   "read_type": "read-perfect",
   "description": "Perfect reads from sample g0_s1 from chrom 1 only"
  },
  "g0": {
   "file_name": "s0.vcf.gz",
   "db": "g0_s0",
   "description": "Sample g0_s0 variants"
  },
  "g1": {
   "file_name": "s1.vcf.gz",
   "graph": "g0_s1",
   "description": "Sample g0_s1 variants"
  },
  "genomedb": {  # Not used by the tools, but used by the analysis
    "file_name": "reddus_genomes.h5",
    "description": "Simulated population"
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
  ("reads", ["s0-r1", "s0-r2", "s0-r3"]),
  ("db", ["g0", "g1"])
  # "db" is not used by bwa but used by bwa-db
  # It is an example of how we can
  # a) add more dimensions to the combinations
  # b) skip a dimension for certain tools
]

benchmark_tools = {
  "tool_analysis": {
    # Description of the tool used to analyze each run. Leave out to indicate no per-tool analysis is done
    # (e.g. for comparative benchmarks where there is no truth set, we go directly to meta-analysis).
    "app_name": "aligner-analysis",
    "tag": "al-anal",
    "inputs": ["bam", "genomedb"],
    # "inputs" will be searched for in outputs of tools and in file_list
    # The name (e.g "bam") must be present in the "output_mapping" field of the tool description
    "outputs": {  # In general, we need to know the suffix for each file we produce.
      "badbam": "bad.bam",
      "perbam": "per.bam",
      "indel-json": "indel.json",
      "mis-plot-circle": "mis.circle.pdf",
      "indel-plot": "indel.pdf"
    }
  },
  "meta_analysis": {
    # Description of the tool used to analyze all runs combined.
    # Leave out to skip this step
    "app_name": "aligner-meta-analysis",
    "tag": "al-meta-anal",
    "inputs": ["indel-json", "mis-plot-circle", "indel-plot", "tool-time"],
    # "inputs" will be searched for in outputs of tools, outputs of tool_analysis and in file_list
    # The name (e.g "tool-time") must be present in the "output_mapping" field of the tool description
    "outputs": {  # In general, we need to know the suffix for each file we produce.
      "result": "zip"
    }
  }
}


# We rename all relevant tool outputs according to our naming convention.
# This mapping indicates what suffix we should use for relevant tool files
# It would be tedious to place this in the output description of each tool
# These are the outputs required and used by the benchmark tools
tool_output_suffix = {
  "bam": "bam",
  "tool-time": "time.txt"
}

bench_spec = bench.create_bench_spec(name='B1', description='Example benchmark',
                                     file_list=file_list, bench_combinations=bench_combinations,
                                     benchmark_tools=benchmark_tools,
                                     tool_output_suffix=tool_output_suffix)


#TODO: How to handle additional files, like indexes, naturally. Ask Boysha
tool_descriptions = [
  ("bwa", {
    "app_name": "bwa-mem",
    "tag": "bwa",
    "inputs": ["ref", "reads"],
    "input_mapping": {
      "ref": "fasta",  # ref is the benchmarking term, fasta is the pin name of the tool input
      "reads": "fastq"
    },
    "output_mapping": {
      "bam": "bwa-bam",  # -> what about index files etc.
      # <benchmarking pin name>: <tool pin name>
      "tool-time": "bwa-time"
    },
    "parameters": {'-p': ''},  # Leave blank to leave defaults
    "description": "Normal BWA"
  }),
  ("poor-bwa", {
    "app_name": "bwa-db",  # This is a bogus BWA pretending to take a VCF also as an input
    "tag": "poor-bwa",
    "input_mapping": {
      "ref": "fasta",
      "reads": "fastq",
      "db": "vcf"
    },
    "output_mapping": {
      "bam": "bwa-db-bam",
      "tool-time": "bwa-db-time"
    },
    "parameters": {'-p': '', '-P': ''},  # Leave blank to leave defaults
    "description": "BWA with no paired end info, pretending to take VCF as input"
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
exe = api.SBGSDK2Executor(bench_run, 'Benchmarking')

#bench.execute_benchmark(bench_run, exe, file_name='bench_state.json', poll_interval=0.1)

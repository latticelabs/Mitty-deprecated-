from mitty.tests import *
from mitty.benchmarking.tool_wrappers import bwa
from mitty.benchmarking import bench
import os
from nose.plugins.skip import SkipTest


def simple_test():
  """BWA bench simple test."""
  try:
    t1 = bwa.Tool(1, 'BWA by Heng Li')
  except RuntimeError, ex:
    raise SkipTest('{:s}. Can not run BWA tool wrapper test'.format(ex))
  file_sets = {1: ('chimera',
                   {'reference_file': os.path.join(example_fasta_genome, 'chimera.fa.gz'),
                    'fastq': null_fastq_name,
                    'fastq_is_paired_interleaved': True})
  }
  tools = {1: ('bwa', t1)}
  tool_p_sets = {1: {1: ('4 threads', {'-t': '4'}), 2: ('8 threads', {'-t': '8'}), 3: ('16 threads', {'-t': '16'})}}
  root_dir = os.path.join(data_dir, 'bench_test')
  bench.run_bench_marks(file_sets, tools, tool_p_sets, [(None, None, 3)],
                        root_dir, bench.analyze_bam)
  assert os.path.exists(os.path.join(root_dir, 'f1_t1_p1/bench.csv'))
  assert os.path.exists(os.path.join(root_dir, 'f1_t1_p1/bench.json'))
  assert os.path.exists(os.path.join(root_dir, 'f1_t1_p2/bench.csv'))
  assert os.path.exists(os.path.join(root_dir, 'f1_t1_p2/bench.json'))
  assert not os.path.exists(os.path.join(root_dir, 'f1_t1_p3/bench.csv'))
  assert not os.path.exists(os.path.join(root_dir, 'f1_t1_p3/bench.json'))
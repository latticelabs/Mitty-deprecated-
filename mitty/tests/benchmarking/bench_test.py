import os

from nose.plugins.skip import SkipTest

from mitty.tests import example_fasta_genome, null_fastq_name, data_dir
from mitty.benchmarking.tool_wrappers import bwa
import mitty.benchmarking.bench as bch


def simple_test():
  """BWA bench simple test."""
  try:
    t1 = bwa.Tool('BWA by Heng Li')
  except RuntimeError, ex:
    raise SkipTest('{:s}. Can not run bench test'.format(ex))
  t1.add_parameter_set('4 threads', {'-t': 4})
  t1.add_parameter_set('8 threads', {'-t': 8})
  t1.add_parameter_set('16 threads', {'-t': 16})

  b = bch.AlignerBenchSets()
  b.add_bench_parameter_set('tight', window=0, block_size=1000000)
  b.add_bench_parameter_set('loose', window=100, block_size=1000000)

  b.add_file_set('chimera', reference_file=os.path.join(example_fasta_genome, 'chimera.fa.gz'),
                 fastq=null_fastq_name, fastq_is_paired_interleaved=True)

  b.add_tool(t1)

  root_dir = os.path.join(data_dir, 'bench_test')
  b.run_bench_marks(root_dir)
  assert os.path.exists(os.path.join(root_dir, 'f0_t0_p0/bps0/bench.csv'))
  assert os.path.exists(os.path.join(root_dir, 'f0_t0_p0/bps1/bench.json'))
  assert os.path.exists(os.path.join(root_dir, 'f0_t0_p1/bps1/bench.csv'))
  assert os.path.exists(os.path.join(root_dir, 'f0_t0_p1/bps0/bench.json'))
  #assert not os.path.exists(os.path.join(root_dir, 'f1_t1_p3'))
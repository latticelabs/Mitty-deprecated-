from mitty.tests import *
from mitty.benchmarking.tool_wrappers import bwa
import os
from nose.plugins.skip import SkipTest


def simple_test():
  """BWA example wrapper simple test."""
  try:
    tool = bwa.Tool(1, 'BWA by Heng Li')
  except RuntimeError, ex:
    raise SkipTest('{:s}. Can not run BWA tool wrapper test'.format(ex))

  file_set = {'reference_file': os.path.join(example_fasta_genome, 'chimera.fa.gz'),
              'fastq': null_fastq_name,
              'fastq_is_paired_interleaved': True}
  tool_p_set = {'-t': '4'}
  tool.setup_files(file_set)
  out_prefix = os.path.join(data_dir, 'bwa_bench_test')
  if not os.path.exists(out_prefix):
    os.makedirs(out_prefix)
  out_file_set = tool.benchmark(file_set, tool_p_set, out_prefix=out_prefix)
  assert os.path.exists(out_file_set['out_bam'])
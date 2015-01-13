import os

from nose.plugins.skip import SkipTest

from mitty.tests import example_fasta_genome, null_fastq_name, data_dir
from mitty.benchmarking.tool_wrappers import bwa


def simple_test():
  """BWA example wrapper simple test."""
  try:
    tool = bwa.Tool('BWA by Heng Li')
  except RuntimeError, ex:
    raise SkipTest('{:s}. Can not run BWA tool wrapper test'.format(ex))

  tool.add_parameter_set('4 threads', {'-t': 4})
  file_set = {'reference_file': os.path.join(example_fasta_genome, 'chimera.fa.gz'),
              'fastq': null_fastq_name,
              'fastq_is_paired_interleaved': True}
  tool.setup_files(file_set)
  out_prefix = os.path.join(data_dir, 'bwa_bench_test')
  if not os.path.exists(out_prefix):
    os.makedirs(out_prefix)
  out_file_set = tool.run(file_set, '4 threads', out_prefix=out_prefix)
  assert os.path.exists(out_file_set['out_bam'])
  #_ = tool.run(file_set, '8 threads', out_prefix=out_prefix)  # This should fail

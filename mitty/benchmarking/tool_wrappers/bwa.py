"""Wrapper for bwa."""
from mitty.benchmarking import bench
import os
import subprocess
import tempfile
import logging
logger = logging.getLogger(__name__)


class Tool(bench.Tool):

  def __init__(self, description='BWA by Heng Li'):
    try:
      subprocess.check_output(['which', 'bwa'])
    except subprocess.CalledProcessError:
      raise RuntimeError('No BWA in path')
    try:
      subprocess.check_output(['which', 'samtools'])
    except subprocess.CalledProcessError:
      raise RuntimeError('No samtools in path')

    super(Tool, self).__init__(name='bwa', description=description)

  def add_parameter_set(self, name, param_dict):
    if name in self.parameter_set:
      raise RuntimeError('This tool parameter set ({:s}) has already been added'.format(str(name)))
    self.parameter_set[name] = param_dict

  def setup_files(self, file_set):
    """The only thing we need to do is make sure that we have indexed the reference."""
    assert os.path.basename(file_set['reference_file']).endswith('fa.gz'), "Reference file must be gzipped fasta (fa.gz)"
    if os.path.exists(file_set['reference_file'] + '.bwt'):
      logger.debug('Bwa index exists')
      return
    logger.debug('Creating bwa index')
    subprocess.call(['bwa', 'index', file_set['reference_file']])

  def run(self, file_set, tool_p_set_name, out_prefix):
    if tool_p_set_name not in self.parameter_set:
      raise RuntimeError('No such tool parameter set ({:s})'.format(tool_p_set_name))
    FNULL = open(os.devnull, 'w')  # Shut up BWA
    temp_sam_fp, temp_sam_fname = tempfile.mkstemp(suffix='.sam')
    # TODO: use the unix time command
    logger.debug('Aligning using bwa. Output to {:s}'.format(temp_sam_fname))
    subprocess.call(['bwa', 'mem', '-v 1'] +
                    ['{:s} {:s}'.format(k, str(v)) for k, v in self.parameter_set[tool_p_set_name].iteritems()] +
                    ['-p' if file_set['fastq_is_paired_interleaved'] else '',
                     file_set['reference_file'], file_set['fastq']], stdout=temp_sam_fp, stderr=FNULL)
    os.close(temp_sam_fp)
    logger.debug('Converting {:s} to bam'.format(temp_sam_fname))
    out_file_set = {'out_bam': os.path.join(out_prefix, 'bwa_bam.bam')}
    with open(out_file_set['out_bam'], 'w') as fp:
      subprocess.call(['samtools', 'view', '-Sb', temp_sam_fname], stdout=fp, stderr=FNULL)
    return out_file_set
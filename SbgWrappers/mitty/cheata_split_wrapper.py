"""This is the wrapper for cheata check

Command line parameters are

cheata check --inbam=INBAM  [-v]
"""
from sbgsdk import define, Process, require
import os


@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class SplitGoodBadAlignments(define.Wrapper):
  class Inputs(define.Inputs):
    bam = define.input(name='BAM', description='aligned .bam file of reads produced by read simulator', required=True)

  class Outputs(define.Outputs):
    good_bam = define.output(name='good BAM', description='.bam (and index) file of correctly aligned reads', list=True)
    bad_bam = define.output(name='bad BAM', description='.bam (and index) file of incorrectly aligned reads', list=True)

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    inbam = self.inputs.bam
    p = Process('python', '/Mitty/cheata.py', 'split', '--inbam', inbam, '-v')
    p.run()

    correct_bam = os.path.splitext(inbam)[0] + '_correct.bam'
    wrong_bam = os.path.splitext(inbam)[0] + '_wrong.bam'

    self.outputs.good_bam.add_file(correct_bam)
    self.outputs.good_bam.add_file(correct_bam + '.bai')

    self.outputs.good_bam.add_file(wrong_bam)
    self.outputs.good_bam.add_file(wrong_bam + '.bai')


def test_checka():
  """Test with the porcine circovirus test data"""

  inputs = {'bam': '/sbgenomics/test-data/sim_reads.bam'}
  params = {}
  wrp = SplitGoodBadAlignments(inputs, params)
  outputs = wrp.test()

  # Test to see the output and indexes exist
  assert os.path.exists(outputs.check_file)

"""This is the wrapper for cheata check

Command line parameters are

cheata check --inbam=INBAM  [-v]
"""
from sbgsdk import define, Process, require
import os


@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class SplitGoodBadAlignments(define.Wrapper):
  class Inputs(define.Inputs):
    bam = define.input(name='BAM', description='aligned .bam file of reads produced by read simulator', required=True)  # , file_types=['.bam'])
    wg = define.input(name='wg.gz', description='Whole genome .gz file')

  class Outputs(define.Outputs):
    good_bam = define.output(name='good BAM', description='.bam (and index) file of correctly aligned reads')  #, list=True)  # , file_types=['.bam', '.bai'])
    bad_bam = define.output(name='bad BAM', description='.bam (and index) file of incorrectly aligned reads')  #, list=True)  # , file_types=['.bam', '.bai'])
    unmapped_bam = define.output(name='unmapped BAM', description='.bam (and index) file of unmapped reads')  #, list=True)  # , file_types=['.bam', '.bai'])

  def execute(self):
    inbam = self.inputs.bam
    wg = self.inputs.wg
    p = Process('python', '/Mitty/cheata.py', 'split', '--inbam', inbam, '--wg', wg, '-v')
    p.run()

    correct_bam = os.path.splitext(inbam)[0] + '_correct.bam'
    bad_bam = os.path.splitext(inbam)[0] + '_wrong.bam'
    unmapped_bam = os.path.splitext(inbam)[0] + '_unmapped.bam'

    self.outputs.good_bam = correct_bam  #.add_file(correct_bam)
    #self.outputs.good_bam.add_file(correct_bam + '.bai')
    self.outputs.good_bam.meta = self.inputs.bam.make_metadata()

    self.outputs.bad_bam = bad_bam  #.add_file(bad_bam)
    #self.outputs.bad_bam.add_file(bad_bam + '.bai')
    self.outputs.bad_bam.meta = self.inputs.bam.make_metadata()

    self.outputs.unmapped_bam = unmapped_bam  # .add_file(unmapped_bam)
    #self.outputs.unmapped_bam.add_file(unmapped_bam + '.bai')
    self.outputs.unmapped_bam.meta = self.inputs.bam.make_metadata()


def test_checka():
  """Test with the porcine circovirus test data"""

  inputs = {'bam': '/sbgenomics/test-data/bwa_aligned_chimera.bam',
            'wg': '/sbgenomics/test-data/chimera.wg.gz'}
  params = {}
  wrp = SplitGoodBadAlignments(inputs, params)
  outputs = wrp.test()

  # Test to see the output and indexes exist
  assert os.path.exists(outputs.good_bam)
  assert os.path.exists(outputs.bad_bam)
  assert os.path.exists(outputs.unmapped_bam)
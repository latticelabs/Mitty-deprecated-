"""This is the wrapper for cheata split

Command line parameters are

cheata split --inbam=INBAM  --wg=WG [-v]
"""
from sbgsdk import define, Process, require
import os


#@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class SplitGoodBadAlignments(define.Wrapper):
  class Inputs(define.Inputs):
    bam = define.input(name='BAM', description='aligned .bam file of reads produced by read simulator', required=True)  # , file_types=['.bam'])
    bai = define.input(name='BAI', description='index for aligned .bam file of reads produced by read simulator', required=True)
    wg = define.input(name='wg.gz', description='Whole genome .gz file', required=True)

  class Outputs(define.Outputs):
    good_bam = define.output(name='good BAM', description='.bam (and index) file of correctly aligned reads')  #, list=True)  # , file_types=['.bam', '.bai'])
    bad_bam = define.output(name='bad BAM', description='.bam (and index) file of incorrectly aligned reads')  #, list=True)  # , file_types=['.bam', '.bai'])
    unmapped_bam = define.output(name='unmapped BAM', description='.bam (and index) file of unmapped reads')  #, list=True)  # , file_types=['.bam', '.bai'])
    data_file = define.output(name='Data file (.pkl)', description='.pkl file of table of read analysis')

  def execute(self):
    # We need the index and the bam file to be in the same directory, and the only way to ensure this is to
    # copy them ourselves
    _inbam = self.inputs.bam
    _inbam_bai = self.inputs.bai
    working_dir = 'WDIR/'
    if not os.path.exists(working_dir):
      os.makedirs(working_dir)
    p = Process('cp', _inbam, working_dir)
    p.run()
    p = Process('cp', _inbam_bai, working_dir)
    p.run()

    inbam = os.path.join(working_dir, os.path.basename(self.inputs.bam))
    wg = self.inputs.wg

    p = Process('python', '/Mitty/cheata.py', 'split', '--inbam', inbam, '--wg', wg, '-v')
    p.run()

    correct_bam = os.path.splitext(inbam)[0] + '_correct.bam'
    bad_bam = os.path.splitext(inbam)[0] + '_wrong.bam'
    unmapped_bam = os.path.splitext(inbam)[0] + '_unmapped.bam'
    data_file = os.path.splitext(inbam)[0] + '_alignment_data.pkl'

    self.outputs.good_bam = correct_bam
    self.outputs.good_bam.meta = self.inputs.bam.make_metadata()

    self.outputs.bad_bam = bad_bam
    self.outputs.bad_bam.meta = self.inputs.bam.make_metadata()

    self.outputs.unmapped_bam = unmapped_bam
    self.outputs.unmapped_bam.meta = self.inputs.bam.make_metadata()

    self.outputs.data_file = data_file


def test_checka():
  """Test with the porcine circovirus test data"""

  inputs = {'bam': '/sbgenomics/test-data/bwa_aligned_chimera.bam',
            'bai': '/sbgenomics/test-data/bwa_aligned_chimera.bam.bai',
            'wg': '/sbgenomics/test-data/chimera.wg.gz'}
  params = {}
  wrp = SplitGoodBadAlignments(inputs, params)
  outputs = wrp.test()

  # Test to see the outputs exist
  assert os.path.exists(outputs.good_bam)
  assert os.path.exists(outputs.bad_bam)
  assert os.path.exists(outputs.unmapped_bam)
  assert os.path.exists(outputs.data_file)
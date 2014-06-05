"""This is the wrapper for cheata check

Command line parameters are

cheata check --inbam=INBAM  --checkfile=CHECK  [-v]
"""
from sbgsdk import define, Process, require
import os


@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class Checka(define.Wrapper):
  class Inputs(define.Inputs):
    bam = define.input(name='BAM file', description='.bam file initially produced by read simulator', required=True)

  class Outputs(define.Outputs):
    check_file = define.output(name='STATS file', description='.pkl file of alignment diagnostics if in "check" mode',
                               required=True)

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    inbam = self.inputs.bam
    check_file = os.path.join(output_dir, os.path.splitext(os.path.basename(inbam))[0] + '_align_diag.pkl')
    p = Process('python', '/Mitty/cheata.py', 'check',
                '--inbam', inbam,
                '--checkfile', check_file,
                '-v')
    p.run()
    self.outputs.check_file = check_file


def test_checka():
  """Test with the porcine circovirus test data"""

  inputs = {'bam': '/sbgenomics/test-data/sim_reads.bam'}
  params = {}
  wrp = Checka(inputs, params)
  outputs = wrp.test()

  # Test to see the output and indexes exist
  assert os.path.exists(outputs.check_file)

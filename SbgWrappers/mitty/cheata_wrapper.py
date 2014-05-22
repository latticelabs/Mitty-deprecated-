"""This is the wrapper for converta

Command line parameters are

cheata --inbam=BAM  --outbam=BAM  --refname=REFNAME  --reflen=REFLEN  [-v]

"""
from sbgsdk import define, Process, require
import os


@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class Cheata(define.Wrapper):
  class Inputs(define.Inputs):
    bam = define.input(name='Simulated BAM', description='.bam file produced by read simulator', required=True)

  class Outputs(define.Outputs):
    bam = define.output(name='Perfectly aligned BAM', description='.bam (and index) file of perfectly aligned reads',
                        list=True)

  class Params(define.Params):
    ref_name = define.string(description='Header for the reference sequence', required=True)
    ref_len = define.integer(description='Length of the reference sequence', required=True)

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    in_file = self.inputs.bam
    out_file = os.path.join(output_dir, os.path.splitext(os.path.basename(in_file))[0] + '_aligned.bam')
    p = Process('python', '/Mitty/cheata.py',
                '--inbam', in_file,
                '--outbam', out_file,
                '--refname', self.params.ref_name,
                '--reflen', self.params.ref_len,
                '-v')
    p.run()
    self.outputs.bam.add_file(out_file)
    self.outputs.bam.add_file(out_file + '.bai')


def test_cheata():
  """Test with the porcine circovirus test data"""

  inputs = {'bam': '/sbgenomics/test-data/sim_reads.bam'}
  params = {'ref_name': 'gi|52547303|gb|AY735451.1|',
            'ref_len': 702}
  wrp = Cheata(inputs, params)
  outputs = wrp.test()

  # Test to see the output and indexes exist
  assert os.path.exists(outputs.bam[0])
  assert os.path.exists(outputs.bam[1])
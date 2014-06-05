"""This is the wrapper for cheata

Command line parameters are

cheata  --inbam=INBAM  --outbam=OUTBAM  --heada=HD  [-v]
"""
from sbgsdk import define, Process, require
import os


@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class Cheata(define.Wrapper):
  class Inputs(define.Inputs):
    bam = define.input(name='Simulated BAM', description='.bam file produced by read simulator', required=True)
    heada = define.input(name='HEADA file', description='.smalla.heada file produced by converta', required=True)

  class Outputs(define.Outputs):
    bam = define.output(name='Perfectly aligned BAM', description='.bam (and index) file of perfectly aligned reads',
                        list=True)

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    inbam = self.inputs.bam
    heada = self.inputs.heada
    outbam = os.path.join(output_dir, os.path.splitext(os.path.basename(inbam))[0] + '_aligned.bam')
    p = Process('python', '/Mitty/cheata.py',
                '--inbam', inbam,
                '--outbam', outbam,
                '--heada', heada,
                '-v')
    p.run()
    self.outputs.bam.add_file(outbam)
    self.outputs.bam.add_file(outbam + '.bai')


def test_cheata():
  """Test with the porcine circovirus test data"""

  inputs = {'bam': '/sbgenomics/test-data/sim_reads.bam',
            'heada': '/sbgenomics/test-data/porcine_circovirus_0.smalla.heada'}
  params = {}
  wrp = Cheata(inputs, params)
  outputs = wrp.test()

  # Test to see the output and indexes exist
  assert os.path.exists(outputs.bam[0])
  assert os.path.exists(outputs.bam[1])
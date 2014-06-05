"""This is the wrapper for converta

Command line parameters are

<fasta> <smalla> [--block_size=BS]

"""
from sbgsdk import define, Process, require
import os


@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class Converta(define.Wrapper):
  class Inputs(define.Inputs):
    fasta = define.input(name='FASTA', description='The input .fasta file', required=True)

  class Outputs(define.Outputs):
    smalla = define.output(name='SMALLA', description='The output .smalla file')
    heada = define.output(name='HEADA', description='The output .smalla.heada file')
    # TODO: when implemented in IGOR, add heada as sidecar file

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    in_file = self.inputs.fasta
    out_file_prefix = os.path.join(output_dir, os.path.splitext(os.path.basename(in_file))[0])
    p = Process('python', '/Mitty/converta.py', in_file, out_file_prefix, '--block_size', 10000000)  # TODO: scale this with allocated resources
    p.run()
    self.outputs.smalla = out_file_prefix  + '.smalla'
    self.outputs.heada = out_file_prefix  + '.smalla.heada'

def test_converta():
  """Test with the porcine circovirus test data"""

  inputs = {'fasta': '/sbgenomics/test-data/porcine_circovirus.fa'}
  params = {}
  wrp = Converta(inputs, params)
  outputs = wrp.test()

  # Test to see the written sequences are correct
  with open(outputs.smalla, 'r') as f1, open('/sbgenomics/test-data/porcine_circovirus.smalla') as f2:
    assert f1.read() == f2.read()
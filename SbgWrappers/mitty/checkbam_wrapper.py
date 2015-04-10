"""This is the wrapper for the reads2bam tool.

Commandline::

  Usage:
    checkbam  --inbam=INBAM  --fout=FOUT  [--block_size=BS] [-v]

  Options:
    --inbam=INBAM           Input bam file name of reads
    --fout=FOUT             Name of output file
    --block_size=BS         Number of reads to process at a time [default: 1000000]
    -v                      Dump detailed logger messages

"""

from sbgsdk import define, Process, require
import json
import os

class Checkbam(define.Wrapper):
  class Inputs(define.Inputs):
    input_bam = define.input(name='Input Bam', required = True)

  class Outputs(define.Outputs):
    output_file = define.output(name='Output file')

  class Params(define.Params):
    output_name = define.string(name='output file name', required=True)
    block_size = define.integer(name='blocksize', default=1000000)

  def execute(self):

    process = Process('perfectbam.py', '--inbam', self.inputs.bam, '--fout', self.params.output_name, '--block_size', self.params.block_size, '-v' )
    process.run()

    # add output files to the output wrapper list
    self.outputs.output_file = self.params.output_name
    self.outputs.output_file.meta = self.inputs.input_bam.make_metadata(file_type='bam')


# test a single snp plugin
def test_bam():
  params = {'bam_name': 'test_bam' }
  inputs = {'genome': ['/Mitty/examples/data/chr1.fa'],
            'input_fastq': '/Mitty/examples/fqs/test_vcf.fq'
           }
  wrp_m = Reads2bam(inputs, params)
  outputs = wrp_m.test()

"""This is the wrapper for the reads2bam tool.

Commandline:

  Usage:
    reads2bam  [-p] --fa_dir=FADIR  --fastq=FASTQ  --bam=BAM  [-v|-V]

  Options:
    -p                      If set indicates fastq is interleaved paired
    --fa_dir=FADIR          Directory where genome is located
    --fasta=FASTQ           Input fasta file
    --bam=BAM               Output bam file
    -v                      Dump detailed logger messages
    -V                      Dump very detailed logger messages


"""

from sbgsdk import define, Process, require
import json
import os
import tarfile

class Reads2bam(define.Wrapper):
  class Inputs(define.Inputs):
    genome = define.input(name='Genome', description='tarball of genome fasta files')
    input_fastq = define.input(name='Input Fastq', required = True)

  class Outputs(define.Outputs):
    output_bam = define.output(name='Bam file', description='A bam file', list=True)

  class Params(define.Params):
    bam_name = define.string(name='output bam name', required=True)
    p = define.boolean(name='Indicates fastq is interleaved pair', default=False)

  def execute(self):

    genome_dir = os.getcwd()
    extract_tar(self.inputs.genome)

    p = '-p' if self.params.p else ''
    process = Process('reads2bam.py', p, '--fa_dir', genome_dir, '--fastq', self.inputs.input_fastq, '--bam', bam_name(self.params.bam_name), '-v' )
    process.run()

    # add output files to the output wrapper list
    self.outputs.output_bam.add_file(bam_name(self.params.bam_name))
    self.outputs.output_bam.add_file(bam_name(self.params.bam_name) + '.bai')
    self.outputs.output_bam.meta = self.inputs.genome.make_metadata(file_type='bam')

def bam_name(name):
  if name.endswith('.bam'):
    return name
  else:
    return name + '.bam'

def extract_tar(tarball):
  tar = tarfile.open(tarball)
  tar.extractall()

# test a single snp plugin
def test_bam():
  params = {'bam_name': 'test_bam' }
  inputs = {'genome': ['/Mitty/examples/data/chr.tar.gz'],
            'input_fastq': '/Mitty/examples/fqs/test_fq.fq'
           }
  wrp_m = Reads2bam(inputs, params)
  outputs = wrp_m.test()

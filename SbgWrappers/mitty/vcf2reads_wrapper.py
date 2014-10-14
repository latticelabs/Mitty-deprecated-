"""This is the wrapper for the vcf2reads tool.

Commandline:

  Usage:
    vcf2reads  --pfile=PFILE  [--corrupt]  [--block_len=BL] [-v|-V]
    vcf2reads plugins
    vcf2reads explain <plugin>

  Options:
    --pfile=PFILE           Name for parameter file
    --corrupt               Write out corrupted reads too.
    --block_len=BL          Consider the sequence in chunks this big. See notes [default: 1000000]
    -v                      Dump detailed logger messages
    -V                      Dump very detailed logger messages
    plugins                 List the available denovo plugins
    explain                 Explain details about the indicated plugin
    <plugin>                The plugin to explain


Parameter file example::

  {
    "files": {
      # An absolute path is left as is
      # a relative path is taken relative to the location of the *script*
      "genome": "/Users/kghose/Data/hg38",   # Directory where genome is located
      "input vcf": "Out/test.vcf",           # Omit or set to none to indicate reads from reference genome
      "output prefix": "Out/reads"           # Output file name prefix
    },
    "rng": {
      "master_seed": 1
    },
    "take reads from": [1,2],           # List the chromosomes the reads should be taken from
    "read_model": "simple_sequential",  # Name of the read plugin to use
    "model_params": {                   # Model specific parameters, need to be under the key "model_params"
      "paired": false,
      "read_len": 100,
      "template_len": 250,
      'read_advance': 20
    }
  }
"""

from sbgsdk import define, Process, require
import json
import os
import tarfile

class Vcf2reads(define.Wrapper):
  class Inputs(define.Inputs):
    genome = define.input(name='Genome', description='tarball of genome fasta files')
    plugins = define.input(name='Plugins', description='The reads plugins')
    input_vcf = define.input(name='Input VCF',
                             description='Omit or set to none to indicate reads from reference genome',
                             default='', list = True)

  class Outputs(define.Outputs):
    output_fastq = define.output(name='Fastq file', description='A fastq file')

  class Params(define.Params):
    fq_name = define.string(name='output fq name', required=True)
    master_seed = define.integer(name='Master seed', default=1)
    chromosomes = define.integer(name='Chromosomes', list = True, description='Chromosomes that reads should be taken from')
    corrupt = define.boolean(name='Corrupt', default=False, description='Write out corrupted reads too')
    block_len = define.integer(name='Block Length', default=1000000, description='Consider the sequence in chunks this big')


  def execute(self):
    genome_dir = os.getcwd()
    extract_tar(self.inputs.genome)
    model_params = json.load(open(self.inputs.plugins, 'r'))
    params_json = {
      "files": {
        "genome": genome_dir,
        "input vcf": find_vcf(self.inputs.input_vcf),
        "output prefix": self.params.fq_name
      },
      "rng": {
        "master_seed": self.params.master_seed
      },
      "take reads from": self.params.chromosomes,
      "read_model": self.inputs.plugins.split('/')[-1].replace('_plugin_params.json',''),
      "model_params": model_params
    }

    corrupt = '--corrupt' if self.params.corrupt else ''
    with open('params.json', 'w') as fp:
      json.dump(params_json, fp, indent=2)
    p = Process('vcf2reads.py', '--pfile', 'params.json', corrupt, '--block_len', self.params.block_len, '-V', )
    p.run()

    # add output files to the output wrapper list
    self.outputs.output_fastq = self.params.fq_name + '.fq'
    self.outputs.output_fastq.meta = self.inputs.genome.make_metadata(file_type='fastq')

def find_vcf(file_list):
  for item in file_list:
    if item.endswith('.vcf.gz'):
      return item

def extract_tar(tarball):
  tar = tarfile.open(tarball)
  tar.extractall()

# test a simple seq plugin
def test_simple_seq():
  from .plugins.reads.simple_seq_wrapper import Simple_Seq

  seq_params = {
  'paired': True,
  'read_len': 100,
  'template_len': 250,
  'read_advance': 25
  }
  seq_inputs = {}
  wrp = Simple_Seq(seq_inputs, seq_params)
  wrp.execute()
  params = {'fq_name': 'test_fq', 'chromosomes': [1] }
  inputs = {'genome': ['/Mitty/examples/data/chr.tar.gz'],
            'plugins': wrp.outputs.json_fragment,
            'input_vcf': ['/Mitty/examples/denovo/Out/test.vcf.gz']
           }
  wrp_m = Vcf2reads(inputs, params)
  outputs = wrp_m.test()

"""This is the wrapper for the denovo tool.

Commandline::

  Usage:
    denovo --pfile=PFILE  [-v]
    denovo plugins
    denovo explain <plugin>

  Options:
    --pfile=PFILE           Name for parameter file
    -v                      Dump detailed logger messages
    plugins                 List the available denovo plugins
    explain                 Explain details about the indicated plugin
    <plugin>                The plugin to explain. If none, explains the parameter file format


Parameter file example::

  {
    "files": {
      "genome": "/Users/kghose/Data/hg38",  # An absolute path is left as is
      "output vcf": "Out/test.vcf"          # a relative path is taken relative to the location of the *script*
    },
    "rng": {
      "master_seed": 1
    },
    "denovo_variant_models": [    # The list of variant models should come under this key
      {
        "snp": {                 # name of the model. To get a list of plugin names type "denovo plugins"
          "chromosome": [1, 2],  # Chromosomes to apply this model to
          "phet": 0.5,           # Parameters required by the model
          "p": 0.01,
          "poisson_rng_seed": 1,
          "base_sub_rng_seed": 2
        }
      },
      {                          # We can chain as many models as we wish
        "delete" : {             # We can repeat models if we want
          "chromosome": [1],
          "model": "delete",
          "phet": 0.5,
          "p_del": 0.01,
          "lam_del": 10
        }
      }
    ]
  }
"""

from sbgsdk import define, Process, require
import json
import os
import tarfile

class Denovo(define.Wrapper):
  class Inputs(define.Inputs):
    genome = define.input(name='Genome', description='tarball of genome fasta files')
    plugins = define.input(name='Plugins', description='The mutation plugins', list=True)

  class Outputs(define.Outputs):
    output_vcf = define.output(name='VCF file', description='A VCF file containing a list of simulated variants', list=True)

  class Params(define.Params):
    vcf_name = define.string(name='output vcf name', required=True)
    master_seed = define.integer(name='Master seed', default=1)
  
  def execute(self):
    genome_dir = os.getcwd()
    extract_tar(self.inputs.genome)
    params_json = {
      "files": {
        "genome": genome_dir,
        "output vcf": self.params.vcf_name
      },
      "rng": {
        "master_seed": self.params.master_seed
      },
      "denovo_variant_models": []
    }
    for in_file in self.inputs.plugins:
      params_json['denovo_variant_models'].append(json.load(open(in_file, 'r')))


    with open('params.json', 'w') as fp:
      json.dump(params_json, fp, indent=2)
    p = Process('denovo.py', '--pfile', 'params.json', '-v', )
    p.run()

    # add output files to the output wrapper list
    self.outputs.output_vcf.add_file(os.path.abspath(self.params.vcf_name+'.vcf'))
    self.outputs.output_vcf.add_file(os.path.abspath(self.params.vcf_name+'.vcf.gz'))
    self.outputs.output_vcf.add_file(os.path.abspath(self.params.vcf_name+'.vcf.gz.tbi'))
    self.outputs.output_vcf.meta = self.inputs.genome.make_metadata(file_type='vcf')

def extract_tar(tarball):
  tar = tarfile.open(tarball)
  tar.extractall()

# test a single snp plugin
def test_simple_snp():
  from .plugins.variants.snp_wrapper import SNP

  snp_params = {
    "model_id": "snp_test",
    "chromosome": [1,2,3,4],      
    "p": 0.001,             
    "phet": 0.5,            
    "base_loc_rng_seed": 1, 
    "base_sub_rng_seed": 2,
    "het_rng_seed": 3,
    "copy_rng_seed": 4
  }

  snp_inputs = {}
  wrp = SNP(snp_inputs, snp_params)
  wrp.execute()

  params = {'vcf_name': 'test_snp' }
  inputs = {'genome': ['/Mitty/examples/data/chr.tar.gz'],
            'plugins': [wrp.outputs.json_fragment]}
  wrp_m = Denovo(inputs, params)
  outputs = wrp_m.test()

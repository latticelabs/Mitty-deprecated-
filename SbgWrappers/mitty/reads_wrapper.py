"""This is the wrapper for reads.py.
Command line parameters are

reads --paramfile=PFILE  [--corrupt]  [--fastq] [--reads_per_block=BL]  [-v]

TODO: We set the flag for --corrupt depending on whether there is a consumer for the corrupt reads.

We set the same number of reads and read range for every sequence in the set.

Example parameter file .json

{
    "whole genome file": "test.wg.gz",
    "whole genome pos file": "test.wg.pos",
    "take reads from": ['1:1', '1:2'],
    "coverage": 5.0,
    "output_file_prefix": "sim_reads",
    "read_model": "tiled_reads",
    "model_params": {
        "paired": false,
        "read_len": 100,
        "template_len": 250,
        "read_advance": 50
    }
}

* If "take reads from" is set to none reads will be generated from all chromosomes. Otherwise reads will be taken only
  from the specified chromosomes
* In the pos file (test.wg.pos in this case) if the pos data for a particular chromosome is not present, we will assume
  that the sequence is the reference sequence
* If the pos file is left empty (None) then we assume this is the reference genome.
"""
from sbgsdk import define, Process, require
import json
import os


#@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class Reads(define.Wrapper):
  class Inputs(define.Inputs):
    wg = define.input(name='wg.gz', description='Whole genome .gz file', required=True)
    wg_pos = define.input(name='.pos.wg.gz', description='Whole genome pos .gz file as produced by vcf2seq')
    plugin = define.input(name='Plugins', description='Hook the output pin of the read plugin you want to use to this.')

  class Outputs(define.Outputs):
    perfect_read_file = define.output(name='Perfect reads',
      description='.bam file containing perfect reads')
    corrupted_read_file = define.output(name='Corrupted reads',
      description='.bam file containing corrupted reads')

  class Params(define.Params):
    coverage = define.real(default=5.0, min=0.0, description='Coverage level', category='General')
    fastq = define.boolean(default=False, description='Use FASTQ as output file format?', category='General')
    #take_reads_from = define.real(default=5.0, min=0.0, description='Coverage level', category='General')

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

    of_prefix = os.path.join(output_dir, os.path.splitext(os.path.basename(self.inputs.wg))[0] + '_sim_reads')

    params_json = {
      "whole genome file": self.inputs.wg,
      "whole genome pos file": self.inputs.wg_pos,
      "take reads from": None,
      "coverage": self.params.coverage,
      "output_file_prefix": of_prefix,
      }
    params_json.update(json.load(open(self.inputs.plugin, 'r')))
    with open('params.json', 'w') as fp:
      json.dump(params_json, fp, indent=2)

    fastq = '--fastq' if self.params.fastq else ''
    p = Process('python', '/Mitty/reads.py', '--paramfile', 'params.json', '--corrupt', fastq, '-v')
    p.run()
    # reads.py produces two files - .bam/.fastq, and _c.bam/_c.fastq

    filetype = 'fastq' if self.params.fastq else 'bam'

    self.outputs.perfect_read_file = params_json['output_file_prefix'] + '.' + filetype
    self.outputs.perfect_read_file.meta = self.outputs.perfect_read_file.make_metadata(file_type=filetype)
    self.outputs.corrupted_read_file = params_json['output_file_prefix'] + '_c.' + filetype
    self.outputs.corrupted_read_file.meta = self.outputs.corrupted_read_file.make_metadata(file_type=filetype)


def test_reads():
  """Test with simple reads from the porcine circovirus test data. Pretend it's a diploid sequence to test input
  of multiple files"""
  json.dump({
    "read_model": "simple_reads",
    "model_params": {
        "paired": True,
        "read_len": 100,
        "template_len": 250,
        "read_loc_rng_seed": 0,
        "error_rng_seed": 1,
        "base_chose_rng_seed": 2,
        "max_p_error": 0.8,
        "k": 0.1
    }
  }, open('/sbgenomics/test-data/read_par.json','w'), indent=2)
  inputs = {'wg': '/sbgenomics/test-data/chimera.wg.gz',
            'plugin': '/sbgenomics/test-data/read_par.json'}
  params = {
    'coverage': 5,
    'fastq': True
  }
  wrp = Reads(inputs, params)
  outputs = wrp.test()

  # we should have two .BAM files
  assert os.path.exists(outputs.perfect_read_file[0])
  assert os.path.exists(outputs.corrupted_read_file[0])
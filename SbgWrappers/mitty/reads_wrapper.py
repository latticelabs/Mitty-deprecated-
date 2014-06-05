"""This is the wrapper for reads.py.
Command line parameters are

[--paramfile=PFILE]  [--corrupt]  [--fastq] [--reads_per_block=BL]  [-v]

TODO: We set the flag for --corrupt depending on whether there is a consumer for the corrupt reads.

We set the same number of reads and read range for every sequence in the set.

Example param file

{
    "input_sequences": ["mutated_1.smalla", "mutated_2.smalla"],
    "total_reads": [100, 100],
    "coverages": [5.0, 5.0]
    "is_this_ref_seq": false,
    "read_ranges": [[0.0, 1.0], [0.0, 1.0]],
    "output_file_prefix": "sim_reads",
    "read_model": "tiled_reads",
    "model_params": {
        "paired": false,
        "read_len": 100,
        "template_len": 250,
        "read_advance": 50
    }
}
"""
from sbgsdk import define, Process, require
import json
import os


@require(mem_mb=4048, cpu=require.CPU_SINGLE)
class Reads(define.Wrapper):
  class Inputs(define.Inputs):
    seq = define.input(name='Input sequence',
                       description='.smalla and .pos file(s) containing the sets of sequences making up the sample.'
                                   ' It is easiest to hook up the output pin of vcf2seq to this pin. Alternatively, '
                                   'pipe any .smalla and .pos files you have to this pin. If this is a reference '
                                   'sequence no need for .pos files.',
                       required=True, list=True)
    plugin = define.input(name='Plugins', description='Hook the output pin of the read plugin you want to use to this.')

  class Outputs(define.Outputs):
    perfect_read_file = define.output(name='Perfect reads',
      description='.bam file containing perfect reads')
    corrupted_read_file = define.output(name='Corrupted reads',
      description='.bam file containing corrupted reads')

  class Params(define.Params):
    coverage = define.real(default=5.0, min=0.0, description='Coverage level', category='General')
    #total_reads = define.integer(default=100, min=1, description='Total number of reads to generate', category='General')
    is_this_ref_seq = define.boolean(default=False,
                                     description='Is this a reference sequence? If true we will not look for .pos file',
                                     category='General')
    read_start = define.real(default=0.0, min=0, max=1, description='From what fraction of the sequence do we start taking reads',
                             category='Advanced')
    read_stop = define.real(default=1.0, min=0, max=1, description='At what fraction of the sequence do we stop taking reads',
                            category='Advanced')

  def execute(self):
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

    #We need to figure out if we are going to have to make up our own output file name
    of_prefix = os.path.splitext(os.path.basename(self.inputs.seq[0]))[0] + '_sim_reads'

    # We only indicate the .smalla files in the json parameters (ignore the .pos files)
    input_smalla = [fname for fname in self.inputs.seq if fname.endswith('.smalla')]
    params_json = {
      "input_sequences": input_smalla,
      #"total_reads": [self.params.total_reads / len(input_smalla)] * len(input_smalla),
      "coverages": [self.params.coverage] * len(input_smalla),
      "is_this_ref_seq": self.params.is_this_ref_seq,
      "read_ranges": [[self.params.read_start, self.params.read_stop]] * len(input_smalla),
      "output_file_prefix": os.path.join(output_dir, of_prefix)
    }
    params_json.update(json.load(open(self.inputs.plugin, 'r')))
    with open('params.json', 'w') as fp:
      json.dump(params_json, fp, indent=2)
    p = Process('python', '/Mitty/reads.py', '--paramfile', 'params.json', '--corrupt')
    p.run()
    # reads.py produces two sets of files - .bam, .bam.bai and _c.bam, _c.bam.bai
    self.outputs.perfect_read_file = params_json['output_file_prefix'] + '.bam'
    self.outputs.perfect_read_file.meta = self.outputs.perfect_read_file.make_metadata(file_type='bam')
    self.outputs.corrupted_read_file = params_json['output_file_prefix'] + '_c.bam'
    self.outputs.corrupted_read_file.meta = self.outputs.corrupted_read_file.make_metadata(file_type='bam')


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
  inputs = {'seq': ['/sbgenomics/test-data/porcine_circovirus_0.smalla', '/sbgenomics/test-data/porcine_circovirus_0.smalla'],
            'plugin': '/sbgenomics/test-data/read_par.json'}
  params = {
    'coverage': 5,
    'is_this_ref_seq': True,
    'read_start': 0.0,
    'read_stop': 1.0
  }
  wrp = Reads(inputs, params)
  outputs = wrp.test()

  # we should have two .BAM files
  assert os.path.exists(outputs.perfect_read_file[0])
  assert os.path.exists(outputs.corrupted_read_file[0])
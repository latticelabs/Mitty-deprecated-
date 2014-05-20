"""This is the wrapper for reads.py.
Command line parameters are

[--paramfile=PFILE]  [--corrupt]  [--fastq] [--reads_per_block=BL]  [-v]

TODO: We set the flag for --corrupt depending on whether there is a consumer for the corrupt reads.

We set the same number of reads and read range for every sequence in the set.

Example param file

{
    "input_sequences": ["mutated_1.smalla", "mutated_2.smalla"],
    "total_reads": [100, 100],
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
                       description='.smalla and .pos file(s) containing the sets of sequences making up the sample',
                       required=True, list=True)
    plugin = define.input(name='Plugins', description='The read plugin')

  class Outputs(define.Outputs):
    perfect_read_file = define.output(name='Perfect reads',
      description='.bam and .bai files containing perfect reads', list=True)
    corrupted_read_file = define.output(name='Corrupted reads',
      description='.bam and .bai files containing corrupted reads', list=True)

  class Params(define.Params):
    total_reads = define.integer(default=100, min=1, description='Number of reads', category='General')
    is_this_ref_seq = define.boolean(default=False,
                                     description='Is this a reference sequence? If true reads will not look for .pos file',
                                     category='General')
    read_start = define.real(default=0.0, min=0, max=1, description='From what fraction of the sequence do we start taking reads',
                             category='Advanced')
    read_stop = define.real(default=1.0, min=0, max=1, description='At what fraction of the sequence do we stop taking reads',
                            category='Advanced')
    output_file_prefix = define.string(default='sim_reads', category='General')

  def execute(self):
    #output_name = self.params.output_vcf_name or input_name + '_variants.vcf'
    #output_name = input_name + '_variants.vcf'
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

    input_smalla = [fname for fname in self.inputs.seq if fname.endswith('.smalla')]  # Discard .pos files
    params_json = {
      "input_sequences": input_smalla,
      "total_reads": [self.params.total_reads] * len(input_smalla),
      "is_this_ref_seq": self.params.is_this_ref_seq,
      "read_ranges": [[self.params.read_start, self.params.read_stop]] * len(input_smalla),
      "output_file_prefix": self.params.output_file_prefix
    }
    params_json.update(json.load(open(self.inputs.plugin, 'r')))
    with open('params.json', 'w') as fp:
      json.dump(params_json, fp, indent=2)
    p = Process('python', '/Mitty/reads.py', '--paramfile', 'params.json', '--corrupt')
    p.run()
    # reads.py produces two files - .bam, and _c.bam
    self.outputs.perfect_read_file = self.params.output_file_prefix + '.bam'
    self.outputs.perfect_read_file.meta = self.outputs.perfect_read_file.make_metadata(file_type='bam')
    self.outputs.corrupted_read_file = self.params.output_file_prefix + '_c.bam'
    self.outputs.corrupted_read_file.meta = self.outputs.corrupted_read_file.make_metadata(file_type='bam')


def test_reads():
  """Test with simple reads from the porcine circovirus test data. Pretend it's a diploid sequence to test input
  of multiple files"""
  json.dump({
    "model": "simple_reads",
    "args": {
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
  inputs = {'seq': ['/sbgenomics/test-data/porcine_circovirus.smalla', '/sbgenomics/test-data/porcine_circovirus.smalla'],
            'plugin': '/sbgenomics/test-data/read_par.json'}
  params = {
    'total_reads': 100,
    'is_this_ref_seq': True,
    'read_start': 0.0,
    'read_stop': 1.0,
    'output_file_prefix': 'sim_reads'
  }
  wrp = Reads(inputs, params)
  outputs = wrp.test()

  # we should have two .BAM files
  assert os.path.exists(outputs.perfect_read_file)
  assert os.path.exists(outputs.corrupted_read_file)
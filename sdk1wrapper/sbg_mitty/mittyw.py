import json
import os

from sbgsdk import define, Process, require
from sbgsdk.file_utils import change_ext


# -------------- Wrappers for Read generator and models ------------------------


@require(mem_mb=8000, cpu=require.CPU_SINGLE)
class Reads(define.Wrapper):
  """This wrapper has to first build a .json file and then run reads on it.

  {
    "files": {
      "reference_file": "reddus_pentalgus.fa.gz",
      "output_prefix": "reads",
      "gzipped": false,
      "dbfile": "reddus_genomes.h5"
    },
    "sample_name" : "g0_s0",
    "rng": {
      "master_seed": 1234567
    },
    "chromosomes": [1, [3, 0.4, 0.6], 4, 5],
    "read_model": "simple_illumina",
    "coverage": 10,
    "coverage_per_block": 1,
    "corrupt": true,
    "variants_only": false,
    "variant_window": 50,
    "model_params": {
      "read_len": 40,
      "template_len_mean": 100,
      "template_len_sd": 20,
      "max_p_error": 0.2,
      "k": 10
    }
  }
  """
  class Inputs(define.Inputs):
    fasta = define.input(required=True, description='Reference FASTA file', name='.fa.gz')
    gdb = define.input(description='Genome HDF5 file', name='.h5')
    read_model = define.input(required=True, description='Read model', name='.json')

  class Outputs(define.Outputs):
    fq_p = define.output(description='Output FASTQ file (perfect)', name='FASTQ (perfect)')
    fq_c = define.output(description='Output FASTQ file (corrupted)', name='FASTQ (corrupted)')

  class Params(define.Params):
    sample_name = define.string(default="g0_s0", description='Name of sample to take read from', required=True)
    rng_master_seed = define.integer(description='Random number seed', min=1, required=True)
    chromosomes = define.integer(description='List of chromosomes to take reads from', list=True, required=True)
    coverage = define.real(default=30, min=0, description='Coverage of reads', required=True)
    corrupt = define.boolean(default=False, description='Generate corrupted reads too?')
    variants_only = define.boolean(default=False, description='Generate reads from neighborhood of variants only?')
    variant_window = define.integer(default=100, min=0, description='Size of neighborhood in bp around variants to take reads from')
    gzipped_fasta = define.boolean(default=False, description='GZIP the fasta files')

  def create_parameter_file(self):
    """Hard coding interleaved to be True for now"""
    read_model_json_fragment = json.load(open(self.inputs.read_model, 'r'))
    read_model_name = read_model_json_fragment.pop('read_model_name_for_sbg_wrappers')
    return {
      "files": {
        "reference_file": self.inputs.fasta,
        "output_prefix": "reads",
        "interleaved": True,
        "gzipped": self.params.gzipped_fasta,
        "dbfile": self.inputs.gdb
      },
      "sample_name": self.params.sample_name,
      "rng": {
        "master_seed": self.params.rng_master_seed
      },
      "chromosomes": self.params.chromosomes,
      "read_model": read_model_name,
      "coverage": self.params.coverage,
      "coverage_per_block": 1,  # Need a good heuristic for this
      "corrupt": self.params.corrupt,
      "variants_only": self.params.variants_only,
      "variant_window": self.params.variant_window,
      "model_params": read_model_json_fragment
    }

  def execute(self):
    with open('empty.fq', 'w') as fp:  # Empty file
      pass
    json.dump(self.create_parameter_file(), open('reads.json', 'w'), indent=2)
    Process('reads', 'generate', '-v', 'reads.json').run()
    self.outputs.fq_p = 'reads.fq.gz' if self.params.gzipped_fasta else 'reads.fq'
    self.outputs.fq_c = ('reads_c.fq.gz' if self.params.gzipped_fasta else 'reads_c.fq') if self.params.corrupt else ('empty.fq')


def test_reads():
  json.dump({"template_len_mean": 250, "max_p_error": 0.01, "k": 20, "read_len": 100, "read_model_name_for_sbg_wrappers": "simple_illumina", "template_len_sd": 30, "gc_bias": {"bias_spread": 0.3, "bias_center": 0.5, "bias_height": 1.5}},
            open('/sbgenomics/test-data/read_model_fragment.json', 'w'))
  inputs = {'fasta': '/sbgenomics/test-data/mock.fa.gz',
            'gdb': '/sbgenomics/test-data/mock_genomes.h5',
            'read_model': '/sbgenomics/test-data/read_model_fragment.json'}
  params = {
    'sample_name': 'g0_s0',
    'rng_master_seed': 10000000,
    'chromosomes': [1, 2, 4],
    'coverage': 1,
    'corrupt': True,
    'variants_only': False,
    'variant_window': 200,
    'gzipped_fasta': False
  }
  outputs = Reads(inputs, params).test()
  assert outputs.fq_p.endswith('reads.fq')
  assert os.path.exists(outputs.fq_p)
  assert os.path.exists(outputs.fq_c)


class SimpleIllumina(define.Wrapper):
  """
  Example parameter set:

  {
    'read_len': 100,          # length of each read
    'template_len_mean': 250, # mean length of template
    'template_len_sd': 30,    # sd of template length
    'max_p_error': 0.01,      # Maximum error rate at tip of read
    'k': 20,                  # exponent
    'gc_bias': {              # Omit if not modeling GC bias
      'bias_center': 0.5,     # Center of uni-modal window
      'bias_height': 1.5,     # Maximum departure from mean coverage (peak of curve)
      'bias_spread': 0.3      # Width measure of curve
    }
  }

  Model defaults: (read_len=100,  template_len_mean=250,  template_len_sd=50,  max_p_error=0.01,  k=20,  gc_bias=None)
  """
  class Outputs(define.Outputs):
    out = define.output(description='JSON fragment for read model', name='ReadModel')

  class Params(define.Params):
    read_len = define.integer(default=100, min=1, description='Length of read')
    template_len_mean = define.integer(default=250, min=1, description='Mean length of template')
    template_len_sd = define.integer(define=30, min=0, description='Sd of template length')
    max_p_error = define.real(default=0.01, description='Maximum error rate at tip of read')
    k = define.real(default=20, description='Exponent of error profile curve')
    has_gc_bias = define.boolean(default=False, description='Set to true if GC bias should be modeled', category='GC bias')
    bias_center = define.real(default=0.5, description='Center of uni-modal window', category='GC Bias')
    bias_height = define.real(default=1.5, description='Maximum departure from mean coverage (peak of curve)', category='GC bias')
    bias_spread = define.real(default=0.3, description='Width measure of curve', category='GC bias')

  def execute(self):
    out_file_name = 'illumina_read_fragment.json'
    json_fragment = {
      'read_model_name_for_sbg_wrappers': 'simple_illumina',
      'read_len': self.params.read_len,
      'template_len_mean': self.params.template_len_mean,
      'template_len_sd': self.params.template_len_sd,
      'max_p_error': self.params.max_p_error,
      'k': self.params.k
    }
    if self.params.has_gc_bias:
      json_fragment['gc_bias'] = {
        'bias_center': self.params.bias_center,
        'bias_height': self.params.bias_height,
        'bias_spread': self.params.bias_spread
      }

    json.dump(json_fragment, open(out_file_name, 'w'))
    self.outputs.out = out_file_name
    #self.outputs.out.meta = self.inputs.fastq.make_metadata()


def test_simple_illumina_read_model():
  """Test with gc bias set"""
  inputs = {}
  params = {
    'read_len': 100,          # length of each read
    'template_len_mean': 250, # mean length of template
    'template_len_sd': 30,    # sd of template length
    'max_p_error': 0.01,      # Maximum error rate at tip of read
    'k': 20,                  # exponent
    'has_gc_bias': True,
    'bias_center': 0.5,     # Center of uni-modal window
    'bias_height': 1.5,     # Maximum departure from mean coverage (peak of curve)
    'bias_spread': 0.3      # Width measure of curve
  }
  outputs = SimpleIllumina(inputs, params).test()
  assert outputs.out.endswith('illumina_read_fragment.json')
  fragment = json.load(open(outputs.out, 'r'))
  assert fragment.get('template_len_sd', None) == 30, fragment
  assert 'gc_bias' in fragment, fragment
  assert fragment['gc_bias'].get('bias_height') == 1.5, fragment

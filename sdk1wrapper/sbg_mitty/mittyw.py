import json
import os
from datetime import datetime
now = datetime.now

from sbgsdk import define, Process, require
from sbgsdk.schema.io_list import IOList
from sbgsdk.file_utils import change_ext

SEED_MAX = (1 << 32) - 1  # Used for seeding rng


def timestamped_prefix(prefix):
  return prefix + '_' + now().strftime('%Y_%m_%d_%H_%M_%S')


# -------------- Wrappers for Genome generator and models ------------------------

# ---- Variant models ----

# ---- delete (mitty.plugins.variants.delete_plugin) ----
class DeletePlugin(define.Wrapper):
  """This is the stock delete plugin. A typical parameter set resembles

    {
      "p": 0.01,           # probability that the deletion will happen at any given base
      "p_end": 0.1,        # probability governing length of deletion
      "min_len": 10,   # Lower bound on deletion lengths
      "max_len": 1000  # upper bound on deletion lengths
    }

    Model defaults: (p=0.01,  p_end=0.1,  min_len=10,  max_len=1000)"""
  class Inputs(define.Inputs):
    dummy = define.input(description='Dummy', required=False)

  class Outputs(define.Outputs):
    out = define.output(description='JSON fragment for delete model', name='VarModel')

  class Params(define.Params):
    p = define.real(default=0.0001, min=0.0, max=1.0, description='Probability that the deletion will happen at any given base')
    p_end = define.real(default=0.1, min=0.0, max=1.0, description='Probability governing length of deletion')
    min_len = define.integer(default=1, min=1, description='Lower bound on deletion lengths')
    max_len = define.integer(default=1, min=1, description='Upper bound on deletion lengths')

  def execute(self):
    out_file_name = 'delete_model_parameter_fragment.json'
    json_fragment = {
      'variant_model_name_for_sbg_wrappers': 'delete',
      "p": self.params.p,           # probability that the deletion will happen at any given base
      "p_end": self.params.p_end,        # probability governing length of deletion
      "min_len": self.params.min_len,   # Lower bound on deletion lengths
      "max_len": self.params.max_len  # upper bound on deletion lengths
    }
    json.dump(json_fragment, open(out_file_name, 'w'))
    self.outputs.out = out_file_name


def test_delete_model():
  """Delete plugin test"""
  inputs = {}
  params = {
    'p': 0.01,
    'p_end': 0.1,
    'min_len': 2,
    'max_len': 10
  }
  outputs = DeletePlugin(inputs, params).test()
  assert outputs.out.endswith('delete_model_parameter_fragment.json')
  fragment = json.load(open(outputs.out, 'r'))
  assert fragment.get('max_len', None) == 10, fragment


# ---- insert (mitty.plugins.variants.insert_plugin) ----
class InsertPlugin(define.Wrapper):
  """Stock insertion model that generates sequences with same base transition matrix as the human genome and creates a
    power-law distribution of insertion lengths.
    A typical parameter set resembles

    {
      "p": 0.0001,        # Per-base probability of having an insertion
      "t_mat": [[ 0.32654629,  0.17292732,  0.24524503,  0.25528135],  # Base transition matrix
                [ 0.3489394,   0.25942695,  0.04942584,  0.3422078],
                [ 0.28778188,  0.21087004,  0.25963262,  0.24171546],
                [ 0.21644706,  0.20588717,  0.24978216,  0.32788362]],
      "p_end": 0.1,       # Probability of chain ending
      "max_len": 1000     # Maximum length of insertion
    }

    Model defaults: (p=0.01,  t_mat=None,  p_end=0.1,  max_len=1000)"""
  class Inputs(define.Inputs):
    dummy = define.input(description='Dummy', required=False)

  class Outputs(define.Outputs):
    out = define.output(description='JSON fragment for insert model', name='VarModel')

  class Params(define.Params):
    p = define.real(default=0.0001, min=0.0, max=1.0, description='Probability that the insertion will happen at any given base')
    p_end = define.real(default=0.1, min=0.0, max=1.0, description='Probability of sequence ending. (Governs length of insertion)')
    t_mat = define.string(default=
                          "[[ 0.32654629,  0.17292732,  0.24524503,  0.25528135],"
                          "[ 0.3489394,   0.25942695,  0.04942584,  0.3422078],"
                          "[ 0.28778188,  0.21087004,  0.25963262,  0.24171546],"
                          "[ 0.21644706,  0.20588717,  0.24978216,  0.32788362]]", description='Base transition matrix')
    max_len = define.integer(default=1, min=1, description='Maximum length of insertion')

  def execute(self):
    out_file_name = 'insert_model_parameter_fragment.json'
    json_fragment = {
      'variant_model_name_for_sbg_wrappers': 'insert',
      "p": self.params.p,
      "p_end": self.params.p_end,
      "t_mat": json.loads(self.params.t_mat),
      "max_len": self.params.max_len
    }
    json.dump(json_fragment, open(out_file_name, 'w'))
    self.outputs.out = out_file_name


def test_insert_model():
  """Insert plugin test"""
  inputs = {}
  params = {
    'p': 0.01,
    'p_end': 0.1,
    't_mat': '[[0.1, 0.2, 0.3, 0.4], [0.12, 0.22, 0.32, 0.42], [0.13, 0.23, 0.33, 0.43], [0.14, 0.24, 0.34, 0.44]]',
    'max_len': 10
  }
  outputs = InsertPlugin(inputs, params).test()
  assert outputs.out.endswith('insert_model_parameter_fragment.json')
  fragment = json.load(open(outputs.out, 'r'))
  assert fragment.get('max_len', None) == 10, fragment
  assert fragment.get('t_mat', None)[1][1] == 0.22, fragment


# ---- snp (mitty.plugins.variants.snp_plugin) ----
class SNPPlugin(define.Wrapper):
  """
  This is the stock SNP plugin. A typical parameter set resembles

  {
    "p": 0.01,              # probability that the SNP will happen at any given base
    "t_mat": [[ 0.32654629,  0.17292732,  0.24524503,  0.25528135],  # Base transition matrix
              [ 0.3489394,   0.25942695,  0.04942584,  0.3422078],   # Leading diagonal is ignored
              [ 0.28778188,  0.21087004,  0.25963262,  0.24171546],
              [ 0.21644706,  0.20588717,  0.24978216,  0.32788362]],
  }

  Model defaults: (p=0.01,  t_mat=None)"""
  class Inputs(define.Inputs):
    dummy = define.input(description='Dummy', required=False)

  class Outputs(define.Outputs):
    out = define.output(description='JSON fragment for SNP model', name='VarModel')

  class Params(define.Params):
    p = define.real(default=0.001, min=0.0, max=1.0, description='Probability that the SNP will happen at any given base')
    t_mat = define.string(default=
                          "[[ 0.32654629,  0.17292732,  0.24524503,  0.25528135],"
                          "[ 0.3489394,   0.25942695,  0.04942584,  0.3422078],"
                          "[ 0.28778188,  0.21087004,  0.25963262,  0.24171546],"
                          "[ 0.21644706,  0.20588717,  0.24978216,  0.32788362]]", description='Base transition matrix')

  def execute(self):
    out_file_name = 'snp_model_parameter_fragment.json'
    json_fragment = {
      'variant_model_name_for_sbg_wrappers': 'snp',
      "p": self.params.p,
      "t_mat": json.loads(self.params.t_mat)
    }
    json.dump(json_fragment, open(out_file_name, 'w'))
    self.outputs.out = out_file_name


def test_snp_model():
  """SNP plugin test"""
  inputs = {}
  params = {
    'p': 0.01,
    't_mat': '[[0.1, 0.2, 0.3, 0.4], [0.12, 0.22, 0.32, 0.42], [0.13, 0.23, 0.33, 0.43], [0.14, 0.24, 0.34, 0.44]]',
  }
  outputs = SNPPlugin(inputs, params).test()
  assert outputs.out.endswith('snp_model_parameter_fragment.json')
  fragment = json.load(open(outputs.out, 'r'))
  assert fragment.get('p', None) == 0.01, fragment
  assert fragment.get('t_mat', None)[2][2] == 0.33, fragment


# ---- uniformdel (mitty.plugins.variants.uniform_deletions) ----
class UniformDeletePlugin(define.Wrapper):
  """
  This deletion model generates uniformly distributed deletion lengths for testing.
  {
    "p": 0.01,           # probability that the deletion will happen at any given base
    "min_len": 10,   # Lower bound on deletion lengths
    "max_len": 1000  # upper bound on deletion lengths
  }

  Model defaults: (p=0.01,  min_len=1,  max_len=1000,  ref=None)"""
  class Inputs(define.Inputs):
    dummy = define.input(description='Dummy', required=False)

  class Outputs(define.Outputs):
    out = define.output(description='JSON fragment for delete model', name='VarModel')

  class Params(define.Params):
    p = define.real(default=0.0001, min=0.0, max=1.0, description='Probability that the deletion will happen at any given base')
    min_len = define.integer(default=1, min=1, description='Lower bound on deletion lengths')
    max_len = define.integer(default=100, min=1, description='Upper bound on deletion lengths')

  def execute(self):
    out_file_name = 'uniform_delete_model_parameter_fragment.json'
    json_fragment = {
      'variant_model_name_for_sbg_wrappers': 'uniformdel',
      "p": self.params.p,           # probability that the deletion will happen at any given base
      "min_len": self.params.min_len,   # Lower bound on deletion lengths
      "max_len": self.params.max_len  # upper bound on deletion lengths
    }
    json.dump(json_fragment, open(out_file_name, 'w'))
    self.outputs.out = out_file_name


def test_uniform_delete_model():
  """Uniform delete plugin test"""
  inputs = {}
  params = {
    'p': 0.01,
    'min_len': 2,
    'max_len': 10
  }
  outputs = UniformDeletePlugin(inputs, params).test()
  assert outputs.out.endswith('uniform_delete_model_parameter_fragment.json')
  fragment = json.load(open(outputs.out, 'r'))
  assert fragment.get('max_len', None) == 10, fragment


# ---- uniformins (mitty.plugins.variants.uniform_insertions) ----
class UniformInsertPlugin(define.Wrapper):
  """
  Insertion model that generates sequences with same base transition matrix as the human genome and creates a
  uniform distribution of insertion lengths. A typical parameter set resembles

  {
    "p": 0.0001,        # Per-base probability of having an insertion
    "t_mat": [[ 0.32654629,  0.17292732,  0.24524503,  0.25528135],  # Base transition matrix
              [ 0.3489394,   0.25942695,  0.04942584,  0.3422078],
              [ 0.28778188,  0.21087004,  0.25963262,  0.24171546],
              [ 0.21644706,  0.20588717,  0.24978216,  0.32788362]],
    "min_len": 2,    # Minimum length of insertion
    "max_len": 30    # Maximum length of insertion
  }

  Model defaults: (p=0.01,  t_mat=None,  min_len=2,  max_len=30)"""
  class Inputs(define.Inputs):
    dummy = define.input(description='Dummy', required=False)

  class Outputs(define.Outputs):
    out = define.output(description='JSON fragment for insert model', name='VarModel')

  class Params(define.Params):
    p = define.real(default=0.0001, min=0.0, max=1.0, description='Probability that the insertion will happen at any given base')
    t_mat = define.string(default=
                          "[[ 0.32654629,  0.17292732,  0.24524503,  0.25528135],"
                          "[ 0.3489394,   0.25942695,  0.04942584,  0.3422078],"
                          "[ 0.28778188,  0.21087004,  0.25963262,  0.24171546],"
                          "[ 0.21644706,  0.20588717,  0.24978216,  0.32788362]]", description='Base transition matrix')
    min_len = define.integer(default=1, min=1, description='Minimum insertion length')
    max_len = define.integer(default=100, min=1, description='Maximum length of insertion')

  def execute(self):
    out_file_name = 'uniform_insert_model_parameter_fragment.json'
    json_fragment = {
      'variant_model_name_for_sbg_wrappers': 'uniformins',
      "p": self.params.p,
      "t_mat": json.loads(self.params.t_mat),
      "min_len": self.params.min_len,
      "max_len": self.params.max_len
    }
    json.dump(json_fragment, open(out_file_name, 'w'))
    self.outputs.out = out_file_name


def test_uniform_insert_model():
  """Insert plugin test"""
  inputs = {}
  params = {
    'p': 0.01,
    'p_end': 0.1,
    't_mat': '[[0.1, 0.2, 0.3, 0.4], [0.12, 0.22, 0.32, 0.42], [0.13, 0.23, 0.33, 0.43], [0.14, 0.24, 0.34, 0.44]]',
    'min_len': 1,
    'max_len': 10
  }
  outputs = UniformInsertPlugin(inputs, params).test()
  assert outputs.out.endswith('uniform_insert_model_parameter_fragment.json')
  fragment = json.load(open(outputs.out, 'r'))
  assert fragment.get('max_len', None) == 10, fragment
  assert fragment.get('t_mat', None)[1][1] == 0.22, fragment


# ---- Site frequency spectrum models ----

# ---- double_exp (mitty.plugins.site_frequency.double_exp) ----
class DoubleExpPlugin(define.Wrapper):
  """
  This stock site frequency plugin employs a double exponential as a simple site frequency model.

  f = a1 * exp(-p * k1) + a2 * exp(-p * k2)   (f is normalized such that sum(f) = 1.0)

  Example parameters:

  {
    "double_exp": {
      "a1": 1.0,        # Amplitude of component 1
      "k1": 20.0,       # decay rate of component 1
      "a2": 0.1,        # Amplitude of component 2
      "k2": 5.0,        # Decay rate of component 2
      "p0": 0.001,      # Lower end of probability range
      "p1": 0.2,        # Upper end of probability range
      "bin_cnt": 30     # Number of bins for probability distribution
    }
  }

  Model defaults: (a1=1.0,  k1=20.0,  a2=0.1,  k2=5.0,  p0=0.001,  p1=0.2,  bin_cnt=21)"""
  class Inputs(define.Inputs):
    dummy = define.input(description='Dummy', required=False)

  class Outputs(define.Outputs):
    out = define.output(description='JSON fragment for site model', name='SiteModel')

  class Params(define.Params):
    a1 = define.real(default=1.0, description='Amplitude of component 1')
    k1 = define.real(default=20.0, description='Decay rate of component 1')
    a2 = define.real(default=1.0, description='Amplitude of component 2')
    k2 = define.real(default=5.0, description='Decay rate of component 2')
    p0 = define.real(default=0.001, description='Lower end of probability range')
    p1 = define.real(default=0.99, description='Upper end of probability range')
    bin_cnt = define.integer(default=30, min=1, description='Bin count for site freq spectrum during re-normalization')

  def execute(self):
    out_file_name = 'double_exp_model_parameter_fragment.json'
    json_fragment = {
      'sfs_model_name_for_sbg_wrappers': 'double_exp',
      "a1": self.params.a1,            # Amplitude of component 1
      "k1": self.params.k1,            # decay rate of component 1
      "a2": self.params.a2,            # Amplitude of component 2
      "k2": self.params.k2,            # Decay rate of component 2
      "p0": self.params.p0,            # Lower end of probability range
      "p1": self.params.p1,            # Upper end of probability range
      "bin_cnt": self.params.bin_cnt   # Number of bins for probability distribution
    }
    json.dump(json_fragment, open(out_file_name, 'w'))
    self.outputs.out = out_file_name


def test_double_exp_model():
  """Double Exp plugin test"""
  inputs = {}
  params = {
      "a1": 1.0,        # Amplitude of component 1
      "k1": 20.0,       # decay rate of component 1
      "a2": 0.1,        # Amplitude of component 2
      "k2": 5.0,        # Decay rate of component 2
      "p0": 0.001,      # Lower end of probability range
      "p1": 0.2,        # Upper end of probability range
      "bin_cnt": 30     # Number of bins for probability distribution
  }
  outputs = DoubleExpPlugin(inputs, params).test()
  assert outputs.out.endswith('double_exp_model_parameter_fragment.json'), outputs.out
  fragment = json.load(open(outputs.out, 'r'))
  assert fragment.get('p1', None) == 0.2, fragment
  assert fragment.get('bin_cnt', None) == 30, fragment


# ---- Population models ----

# ---- standard (mitty.plugins.population.standard) ----
class StandardPopPlugin(define.Wrapper):
  """
  Standard population model that picks variants randomly from the master list based on probability value
  Example parameters:

  {
    "standard": {
      "sample_size": 10,
      "force_homozygous": True, # Force sample variants to be 1|1
      "filter_multi_allele": False,  # Take out locii with different variants on the two copies
      "filter_hom": False,  # Take out homozygous
      "max_v_count": 100,  # Maximum number of variants
      "min_v_spacing": 1000  # Minimum gap between variants on same copy
    }
  }

  Model defaults: (sample_size=10,  force_homozygous=False,  filter_multi_allele=False,  filter_hom=False,  max_v_count=None,  min_v_spacing=None)"""
  class Inputs(define.Inputs):
    dummy = define.input(description='Dummy', required=False)

  class Outputs(define.Outputs):
    out = define.output(description='JSON fragment for population model', name='PopModel')

  class Params(define.Params):
    sample_size = define.integer(default=10, min=1, description='Size of sample')
    force_homozygous = define.boolean(default=False, description='Force sample variants to be 1|1')
    filter_multi_allele = define.boolean(default=False, description='Take out locii with different variants on the two copies')
    filter_hom = define.boolean(default=False, description='Take out homozygous')
    max_v_count = define.integer(default=None, description='Maximum number of variants')
    min_v_spacing = define.integer(default=None, description='Minimum gap between variants on same copy')

  def execute(self):
    out_file_name = 'standard_pop_model_parameter_fragment.json'
    json_fragment = {
      'population_model_name_for_sbg_wrappers': 'standard',
      'sample_size': self.params.sample_size,
      'force_homozygous': self.params.force_homozygous, # Force sample variants to be 1|1
      'filter_multi_allele': self.params.filter_multi_allele,  # Take out locii with different variants on the two copies
      'filter_hom': self.params.filter_hom  # Take out homozygous
    }
    if self.params.max_v_count:
      json_fragment['max_v_count'] = self.params.max_v_count  # Maximum number of variants
    if self.params.min_v_spacing:
      json_fragment['min_v_spacing'] = self.params.min_v_spacing  # Minimum gap between variants on same copy

    json.dump(json_fragment, open(out_file_name, 'w'))
    self.outputs.out = out_file_name


def test_standard_pop_model():
  """Standard pop model plugin test"""
  inputs = {}
  params = {
    "sample_size": 10,
    "force_homozygous": True,  # Force sample variants to be 1|1
    "filter_multi_allele": False,  # Take out locii with different variants on the two copies
    "filter_hom": False,  # Take out homozygous
    "max_v_count": 100,  # Maximum number of variants
    "min_v_spacing": None  # Minimum gap between variants on same copy
  }
  outputs = StandardPopPlugin(inputs, params).test()
  assert outputs.out.endswith('standard_pop_model_parameter_fragment.json')
  fragment = json.load(open(outputs.out, 'r'))
  assert fragment.get('sample_size', None) == 10, fragment
  assert fragment.get('force_homozygous', None) is True, fragment
  assert fragment.get('filter_hom', None) is False, fragment
  assert 'min_v_spacing' not in fragment, fragment


# ---- vn (mitty.plugins.population.vn) ----
class VnPopPlugin(define.Wrapper):
  """
  A population model that creates samples with more and more variants. Suitable for the aligner paper experiments
  ^ = intersection
  E = subset

  vx ^ v0 = v0
  vx ^ v1 = v0
  ...
  vx ^ vn = v0

  v0 E v1
  v1 E v2
  v2 E v3
  ...
  v(n-1) E vn

  This plugin does not honor the site frequency spectrum model and ignores the original 'p' values

  {
    "standard": {
      "p_vx": 0.2,
      "p_vn": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
    }
  }"""
  class Inputs(define.Inputs):
    dummy = define.input(description='Dummy', required=False)

  class Outputs(define.Outputs):
    out = define.output(description='JSON fragment for population model', name='PopModel')

  class Params(define.Params):
    p_vx = define.real(min=0.0, max=1.0, description='p_vx')
    p_vn = define.real(min=0.0, max=1.0, list=True, description='p_vn')

  def execute(self):
    out_file_name = 'vn_pop_model_parameter_fragment.json'
    json_fragment = {
      'population_model_name_for_sbg_wrappers': 'vn',
      'p_vx': self.params.p_vx,
      'p_vn': self.params.p_vn
    }
    json.dump(json_fragment, open(out_file_name, 'w'))
    self.outputs.out = out_file_name


def test_vn_pop_model():
  """Vn pop model plugin test"""
  inputs = {}
  params = {
      "p_vx": 0.2,
      "p_vn": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
  }
  outputs = VnPopPlugin(inputs, params).test()
  assert outputs.out.endswith('vn_pop_model_parameter_fragment.json'), outputs.out
  fragment = json.load(open(outputs.out, 'r'))
  assert fragment.get('p_vx', None) == 0.2, fragment
  assert fragment.get('p_vn', None) == [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], fragment


# -------------- Wrapper for Genome generator ----------------------------------
@require(mem_mb=8000, cpu=require.CPU_SINGLE)
class Genomes(define.Wrapper):
  """The wrapper takes variant models, site spectrum models and population models, builds a parameter file and then
  executes based on that parameter file

  {
    # Path notes: an absolute path is left as is. A relative path is taken relative to the parameter file location
    "files": {
      "reference_file": "/Users/kghose/Data/hg38/hg38.fa.gz",
      "dbfile": "Out/test.db"  # Output database file
    },
    "rng": {
      "master_seed": 1
    },
    "population_model": {
      "standard": {            # Name of sample chooser (population) model.
        "sample_size": 1,      # population model parameters
      }
    },
    "site_model": {
        "double_exp": {    # Name of model that handles the site frequency spectrum
          "k1": 0.1,       # spectrum model parameters
          "k2": 2.0,
          "p0": 0.001,
          "p1": 0.2,
          "bin_cnt": 30
        }
      }
    "chromosomes": [1, 2]  # Chromosomes to apply the models to
    "variant_models": [    # The list of variant models should come under this key
      {
        "snp": {           # name of the model.
          "p": 0.01        # Parameters required by the model
        }
      },
      {                    # We can chain as many models as we wish. We can also repeat models.
        "delete" : {
          "p": 0.01
        }
      }
    ]
  }

  """
  class Inputs(define.Inputs):
    fasta = define.input(description='Reference FASTA file', name='Ref', required=True)
    population_model = define.input(description='Population model', name='PopModel', required=True)
    site_model = define.input(description='Site model', name='SiteModel')
    variant_models = define.input(description='Variant models', name='VarModels', list=True, required=True)

  class Outputs(define.Outputs):
    gdb = define.output(description='Genome HDF5 file', name='GenomeDB')
    vcf_files = define.output(description='Selected VCF files from database', name='VCFs', list=True)

  class Params(define.Params):
    rng_master_seed = define.integer(description='Random number seed', min=1, max=SEED_MAX, required=True)
    prefix = define.string(default='mitty', description='Prefix to add to data file(s) (date/time stamp will be added automatically)')
    chromosomes = define.integer(description='List of chromosomes to generate variants on', list=True, required=True)
    sample_names = define.string(description='Name of samples to write out VCFs for', list=True)

  @staticmethod
  def parse_json_fragment(wrapper_code_key, json_fragment_fname):
    """Given a model json fragment, extract the required parameter snippet"""
    json_fragment = json.load(open(json_fragment_fname, 'r'))
    if wrapper_code_key in json_fragment:
      model_name = json_fragment.pop(wrapper_code_key)
      this_fragment = {model_name: json_fragment}
    else:
      raise RuntimeError('Unexpected key in json fragment ' + str(json_fragment))
    return this_fragment

  @staticmethod
  def parse_model(wrapper_code_key, model_category, json_fragment_fname):
    """Given a model json fragment, extract the required parameter snippet"""
    return {
      model_category: Genomes.parse_json_fragment(wrapper_code_key, json_fragment_fname) if type(json_fragment_fname) is not IOList
      else [Genomes.parse_json_fragment(wrapper_code_key, jf) for jf in json_fragment_fname]
    }

  def create_parameter_file(self, prefix):
    """Given the inputs, parameters and models, create a parameter file"""
    params = {
      "files": {
        "reference_file": self.inputs.fasta,
        "dbfile": prefix + "_genomes.h5"
      },
      "rng": {
        "master_seed": self.params.rng_master_seed
      },
      "chromosomes": self.params.chromosomes,
    }
    params.update(Genomes.parse_model('population_model_name_for_sbg_wrappers',
                                      'population_model', self.inputs.population_model))  # population_model
    params.update(Genomes.parse_model('sfs_model_name_for_sbg_wrappers',
                                      'site_model', self.inputs.site_model))  # site_model
    params.update(Genomes.parse_model('variant_model_name_for_sbg_wrappers',
                                      'variant_models', self.inputs.variant_models))  # variant_models
    return params

  def execute(self):
    prefix = timestamped_prefix(self.params.prefix)
    json.dump(self.create_parameter_file(prefix), open('variants.json', 'w'), indent=2)
    Process('genomes', 'generate', '-v', 'variants.json').run()
    self.outputs.gdb = prefix + '_genomes.h5'
    for n, sample in enumerate(self.params.sample_names):
      Process('genomes', 'genome-file', 'write-vcf', '--sample-name', sample, prefix + '_genomes.h5', prefix + '_{:d}.vcf'.format(n)).run()
    self.outputs.vcf_files = [prefix + '_{:d}.vcf'.format(n) for n in range(len(self.params.sample_names))]


def run_standard_pop_model():
  """Standard pop model plugin test"""
  inputs = {}
  params = {
    "sample_size": 10,
    "force_homozygous": True,  # Force sample variants to be 1|1
    "filter_multi_allele": False,  # Take out locii with different variants on the two copies
    "filter_hom": False,  # Take out homozygous
    "max_v_count": 100,  # Maximum number of variants
    "min_v_spacing": None  # Minimum gap between variants on same copy
  }
  return StandardPopPlugin(inputs, params).test()


def run_double_exp_model():
  """Run double exp plugin"""
  inputs = {}
  params = {
      "a1": 1.0,        # Amplitude of component 1
      "k1": 20.0,       # decay rate of component 1
      "a2": 0.1,        # Amplitude of component 2
      "k2": 5.0,        # Decay rate of component 2
      "p0": 0.001,      # Lower end of probability range
      "p1": 0.2,        # Upper end of probability range
      "bin_cnt": 30     # Number of bins for probability distribution
  }
  return DoubleExpPlugin(inputs, params).test()


def run_snp_model():
  inputs = {}
  params = {
    'p': 0.01,
    't_mat': '[[0.1, 0.2, 0.3, 0.4], [0.12, 0.22, 0.32, 0.42], [0.13, 0.23, 0.33, 0.43], [0.14, 0.24, 0.34, 0.44]]',
  }
  return SNPPlugin(inputs, params).test()


def run_insert_model():
  inputs = {}
  params = {
    'p': 0.01,
    'p_end': 0.1,
    't_mat': '[[0.1, 0.2, 0.3, 0.4], [0.12, 0.22, 0.32, 0.42], [0.13, 0.23, 0.33, 0.43], [0.14, 0.24, 0.34, 0.44]]',
    'max_len': 10
  }
  return InsertPlugin(inputs, params).test()


def test_genomes():
  """Test `genomes` wrapper"""
  population_model = run_standard_pop_model()
  site_model = run_double_exp_model()
  v1_model, v2_model = run_snp_model(), run_insert_model()

  inputs = {
    'fasta': '/sbgenomics/test-data/mock.fa.gz',
    'population_model': population_model.out,
    'site_model': site_model.out,
    'variant_models': [v1_model.out, v2_model.out]
  }
  params = {
    'rng_master_seed': 10000000,
    'chromosomes': [1, 2, 4],
    'sample_names': ['g0_s0', 'g0_s1', 'g0_s2']
  }
  outputs = Genomes(inputs, params).test()
  assert outputs.gdb.endswith('genomes.h5')
  assert len(outputs.vcf_files) == 3


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
    fasta = define.input(required=True, description='Reference FASTA file', name='Ref')
    gdb = define.input(description='Genome HDF5 file', name='GenomeDB')
    read_model = define.input(required=True, description='Read model', name='ReadModel')

  class Outputs(define.Outputs):
    fq_p = define.output(description='Output FASTQ file (perfect)', name='FASTQ (perfect)')
    fq_c = define.output(description='Output FASTQ file (corrupted)', name='FASTQ (corrupted)')

  class Params(define.Params):
    sample_name = define.string(default="g0_s0", description='Name of sample to take read from', required=True)
    rng_master_seed = define.integer(description='Random number seed', min=1, required=True)
    prefix = define.string(default='mitty', description='Prefix to add to data file(s) (date/time stamp will be added automatically)')
    chromosomes = define.integer(description='List of chromosomes to take reads from', list=True, required=True)
    coverage = define.real(default=30, min=0, description='Coverage of reads', required=True)
    corrupt = define.boolean(default=False, description='Generate corrupted reads too?')
    variants_only = define.boolean(default=False, description='Generate reads from neighborhood of variants only?')
    variant_window = define.integer(default=100, min=0, description='Size of neighborhood in bp around variants to take reads from')
    gzipped_fasta = define.boolean(default=False, description='GZIP the fasta files')

  def create_parameter_file(self, prefix):
    """Hard coding interleaved to be True for now"""
    read_model_json_fragment = json.load(open(self.inputs.read_model, 'r'))
    read_model_name = read_model_json_fragment.pop('read_model_name_for_sbg_wrappers')
    return {
      "files": {
        "reference_file": self.inputs.fasta,
        "output_prefix": prefix + "_reads",
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
    prefix = timestamped_prefix(self.params.prefix)
    json.dump(self.create_parameter_file(prefix), open('reads.json', 'w'), indent=2)
    Process('reads', 'generate', '-v', 'reads.json').run()
    self.outputs.fq_p = (prefix + '_reads.fq.gz') if self.params.gzipped_fasta else (prefix + '_reads.fq')
    self.outputs.fq_c = ((prefix + '_reads_c.fq.gz') if self.params.gzipped_fasta else (prefix + '_reads_c.fq')) if self.params.corrupt else None


def test_reads():
  """Perfect and corrupted reads"""
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


def test_reads2():
  """Only perfect reads"""
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
    'corrupt': False,
    'variants_only': False,
    'variant_window': 200,
    'gzipped_fasta': False
  }
  outputs = Reads(inputs, params).test()
  assert outputs.fq_p.endswith('reads.fq')
  assert os.path.exists(outputs.fq_p)


# ---- Read models -------

# -- Simple Illumina plugin --
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
  class Inputs(define.Inputs):
    dummy = define.input(description='Dummy', required=False)

  class Outputs(define.Outputs):
    out = define.output(description='JSON fragment for read model', name='ReadModel')

  class Params(define.Params):
    read_len = define.integer(default=100, min=1, description='Length of read')
    template_len_mean = define.integer(default=250, min=1, description='Mean length of template')
    template_len_sd = define.integer(default=20, min=0, description='Sd of template length')
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


# --------- Benchmarking tools wrappers ------------
class Perfectbam(define.Wrapper):
  class Inputs(define.Inputs):
    inbam = define.input(required=True, description='BAM file to analyze', name='BAM')

  class Outputs(define.Outputs):
    bad_bam = define.output(description='BAM with misaligned reads', name='BADBAM')
    per_bam = define.output(description='BAM with all reads correctly aligned and alignment analysis results set', name='PERBAM')

  class Params(define.Params):
    window = define.integer(default=100, min=0, description='Size of tolerance window for marking alignments correct')
    perfect_bam = define.boolean(default=False, description='Perfect BAM has full read information (results in large size)')
    cigar_errors = define.boolean(default=False, description='CIGAR errors result in reads being classified as misaligned')
    x = define.boolean(default=False, description='Use extended CIGAR ("X"s and "="s) rather than traditional CIGAR (just "M"s)')

  def execute(self):
    bad_bam, per_bam = change_ext(self.inputs.inbam, 'bad.bam'), change_ext(self.inputs.inbam, 'per.bam')
    Process('samtools', 'sort', '-@', 8, '-f', self.inputs.inbam, 'sorted.bam').run()
    Process('samtools', 'index', 'sorted.bam').run()
    argument_list = ['perfectbam', '-v', '--window', self.params.window] + \
                    (['--perfect-bam'] if self.params.perfect_bam else []) + \
                    (['--cigar-errors'] if self.params.cigar_errors else []) + \
                    (['-x'] if self.params.x else []) + \
                    ['sorted.bam', '--per-bam', per_bam, '--bad-bam', bad_bam, '--no-index']
    Process(*argument_list).run()
    self.outputs.bad_bam = bad_bam
    self.outputs.per_bam = per_bam


def test_perfectbam():
  """Test perfectbam"""
  inputs = {
    'inbam': '/sbgenomics/test-data/mock.bam'
  }
  params = {
    'window': 100,
    'perfect_bam': False,
    'cigar_errors': False,
    'x': False
  }
  outputs = Perfectbam(inputs, params).test()
  assert os.path.exists(outputs.bad_bam)
  assert os.path.exists(outputs.per_bam)


@require(cpu=require.CPU_ALL)
class Alindel(define.Wrapper):
  """Wraps BAM sort/index and alindel"""
  class Inputs(define.Inputs):
    per_bam = define.input(required=True, description='Perfect BAM file', name='PERBAM')
    gdb = define.input(description='Genome HDF5 file', name='.h5')

  class Outputs(define.Outputs):
    out_json = define.output(description='.json file with indel accuracy analysis', name='IndelAnal')

  class Params(define.Params):
    indel_range = define.integer(default=100, min=0, description='Indel range to compute')
    sample_name = define.string(description='Sample name')

  def execute(self):
    indel_json = change_ext(self.inputs.per_bam, 'indel.json')
    # First we have to sort and index the BAM
    Process('samtools', 'sort', '-@', 8, '-f', self.inputs.per_bam, 'sorted.bam').run()
    Process('samtools', 'index', 'sorted.bam').run()
    argument_list = ['alindel', '--indel-range', self.params.indel_range] + \
                    (['--sample-name', self.params.sample_name] if self.params.sample_name else []) + \
                    ['sorted.bam', self.inputs.gdb, indel_json]
    Process(*argument_list).run()
    self.outputs.out_json = indel_json


def test_alindel():
  """Test alindel"""
  inputs = {
    'per_bam': '/sbgenomics/test-data/mock_per.bam',
    'gdb': '/sbgenomics/test-data/mock_genomes.h5'
  }
  params = {
    'indel_range': 20,
    'sample_name': 'g0_s0'
  }
  outputs = Alindel(inputs, params).test()
  assert os.path.exists(outputs.out_json)
  data = json.load(open(outputs.out_json, 'r'))
  assert 'templates_within_feature_but_read_outside' in data, data


class AlindelPlot(define.Wrapper):
  """  -f PATH                Indel file(s) to plot
  -o PATH                Output file name
  -l TEXT                Label to go with each file
  --win INTEGER          Size of median filter window to smooth plots
  --indel-range INTEGER  Maximum indel length to plot
  --title TEXT           Title
  --help                 Show this message and exit."""
  class Inputs(define.Inputs):
    indel_json = define.input(required=True, list=True, description='List of indel.json files to plot', name='IndelAnal')

  class Outputs(define.Outputs):
    figure_file = define.output(description='Output figure', name='IndelPlot')

  class Params(define.Params):
    #TODO: Smart labels?
    prefix = define.string(default='mitty', description='Prefix to add to data file(s) (date/time stamp will be added automatically)')
    window = define.integer(default=5, min=0, description='Smoothing window to apply to plot')
    indel_range = define.integer(default=100, min=0, description='Indel range to show')
    plot_title = define.string(default='Alignment accuracy', description='Plot title')
    pdf_plot = define.boolean(default=True, description='Plot as pdf or png')

  def execute(self):
    prefix = timestamped_prefix(self.params.prefix)
    out_plot_name = prefix + '_indel_plot.pdf' if self.params.pdf_plot else prefix + '_indel_plot.png'
    argument_list = ['alindel_plot', '-o', out_plot_name, '--win', self.params.window,
                     '--indel-range', self.params.indel_range, '--title', self.params.plot_title]
    for fname in self.inputs.indel_json:
      argument_list += ['-f', fname]
    Process(*argument_list).run()
    self.outputs.figure_file = out_plot_name


def test_alindel_plot():
  """Test alindel plot"""
  inputs = {
    'indel_json': ['/sbgenomics/test-data/mock_indel.json', '/sbgenomics/test-data/mock_indel.json']
  }
  params = {
    'window': 3,
    'indel_range': 0,
    'plot_title': 'All your gnome belong to us',
    'pdf_plot': True
  }
  outputs = AlindelPlot(inputs, params).test()
  assert os.path.exists(outputs.figure_file)


class MisalignmentPlot(define.Wrapper):
  """
      --circle PATH           Name of figure file for circle plot
      --matrix PATH           Name of figure file for matrix plot
      --bin-size FLOAT        Bin size in Mb
      --scaling-factor FLOAT  Scale size of disks/lines in plot"""
  class Inputs(define.Inputs):
    bad_bam = define.input(required=True, description='BAD BAM produced by perfectbam', name='BADBAM')

  class Outputs(define.Outputs):
    circle_figure_file = define.output(description='Circle figure plot', name='CirclePlot')
    matrix_figure_file = define.output(description='Matrix figure plot', name='MatrixPlot')

  class Params(define.Params):
    bin_size = define.real(default=1, min=0, description='Bin size in Mb')
    scaling_factor = define.real(default=1, min=0, description='Scale size of disks/lines in plot')
    pdf_plot = define.boolean(default=True, description='Plot as pdf or png')

  def execute(self):
    # First we have to sort and index the BAM
    name = os.path.basename(self.inputs.bad_bam)
    base, old_ext = os.path.splitext(name)
    ext = 'pdf' if self.params.pdf_plot else 'png'
    circle_plot_name = base + '_circle.' + ext
    matrix_plot_name = base + '_matrix.' + ext

    # Sort and index input file
    Process('samtools', 'sort', '-@', 8, '-f', self.inputs.bad_bam, 'sorted.bam').run()
    Process('samtools', 'index', 'sorted.bam').run()
    argument_list = ['misplot', '--circle', circle_plot_name, '--matrix', matrix_plot_name,
                     '--bin-size', self.params.bin_size, '--scaling-factor', self.params.scaling_factor,
                     'sorted.bam']
    Process(*argument_list).run()
    self.outputs.circle_figure_file = circle_plot_name
    self.outputs.matrix_figure_file = matrix_plot_name


def test_misalignment_plot():
  """Test mis-alignment plot"""
  inputs = {
    'bad_bam': '/sbgenomics/test-data/mock_bad.bam'
  }
  params = {
    'bin_size': 0.01,
    'scaling_factor': 1.0,
    'pdf_plot': True
  }
  outputs = MisalignmentPlot(inputs, params).test()
  assert os.path.exists(outputs.circle_figure_file)
  assert os.path.exists(outputs.matrix_figure_file)
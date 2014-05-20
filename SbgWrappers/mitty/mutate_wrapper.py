"""This is the wrapper for the mutate tool.
Command line parameters are

--paramfile=PFILE  [--block_size=BS] [-v]

Input json file for mutate.py is like

{
    "input smalla file": "Data/porcine_circovirus.smalla",
    "output vcf file": "TEST-DATA/variants.vcf",
    "chromosome": "1",
    "mutations": {
        "snp": {
            "model": "snp",
            "start_snps_frac": 0.1,
            "stop_snps_frac":  0.3,
            "phet": 0,
            "p": 0.01,
            "poisson_rng_seed": 1,
            "base_sub_rng_seed": 2
        },
        "delete": {
            "model": "delete",
            "start_dels_frac": 0.4,
            "stop_dels_frac":  0.6,
            "phet": 0,
            "p_del": 0.01,
            "lam_del": 10,
            "del_loc_rng_seed": 10,
            "del_len_rng_seed": 1
        },
        "insert": {
            "model": "insert",
            "start_ins_frac": 0.7,
            "stop_ins_frac":  0.9,
            "phet": 0,
            "p_ins": 0.01,
            "lam_ins": 10,
            "ins_loc_rng_seed": 0,
            "ins_len_rng_seed": 1,
            "base_sel_rng_seed": 2
        }
    }
}

We get pieces of the .json file corresponding to entries under "mutations" from the plugin wrappers and then piece it
together here.
"""
from sbgsdk import define, Process, require
import json
import os


@require(mem_mb=1024, cpu=require.CPU_SINGLE)
class Mutate(define.Wrapper):
  class Inputs(define.Inputs):
    ref = define.input(name='Reference', description='The reference sequence (.smalla format)')
    plugins = define.input(name='Plugins', description='The mutation plugins', list=True)

  class Outputs(define.Outputs):
    vcf = define.output(name='VCF file', description='A VCF file containing a list of simulated variants', list=True)

  class Params(define.Params):
    vcf_file_name = define.string(default='', description='Output file name. Leaving this blank will generate a file name with _variants.vcf added to the sequence name')
    chromosome = define.string(default='1', description='Chromosome number')

  def execute(self):
    input_name = os.path.splitext(os.path.basename(self.inputs.ref))[0]
    output_basename = input_name + '_variants.vcf' if self.params.vcf_file_name == '' else self.params.vcf_file_name
    output_dir = 'OUTPUT'
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    output_absolute_path = os.path.join(output_dir, output_basename)
    mutations = {}
    for in_file in self.inputs.plugins:
      mutations.update(json.load(open(in_file, 'r')))

    params_json = {
      "input smalla file": self.inputs.ref,
      "output vcf file": output_absolute_path,
      "chromosome": self.params.chromosome,
      "mutations": mutations
    }
    with open('params.json', 'w') as fp:
      json.dump(params_json, fp, indent=2)
    p = Process('python', '/Mitty/mutate.py', '--paramfile', 'params.json', '-v')
    p.run()
    # mutate.py produces three files - the .vcf, the gzipped form .vcf.gz and the tabix index .vcf.gz.tbi
    self.outputs.vcf.add_file(output_absolute_path)
    self.outputs.vcf.add_file(output_absolute_path + '.gz')
    self.outputs.vcf.add_file(output_absolute_path + '.gz.tbi')
    self.outputs.vcf.meta = self.inputs.ref.make_metadata(file_type='vcf')


def test_mutate_simple():
  """Test with one plugin."""
  from .plugins.mutation.snp_wrapper import SNP  # Thanks Nebojsa

  snp_params = {
          "model_id": "snp_test",
          "start_snps_frac": 0.1,
          "stop_snps_frac": 0.3,
          "phet": 0,
          "p": 0.01,
          "het_rng_seed": 3,
          "strand_rng_seed": 4,
          "poisson_rng_seed": 1,
          "base_sub_rng_seed": 2
  }
  snp_inputs = {}
  wrp = SNP(snp_inputs, snp_params)
  wrp.execute()
  params = {'chromosome': '1'}
  inputs = {'ref': '/sbgenomics/test-data/porcine_circovirus.smalla',
            'plugins': [wrp.outputs.json_fragment]}
  wrp_m = Mutate(inputs, params)
  outputs = wrp_m.test()

  expected_line = '1	169	.	A	T	96	PASS	.	GT	1/1'
  with open(outputs.vcf[0], 'r') as f:
    assert f.readlines()[5].strip() == expected_line


def test_mutate_plugin_list_model_ids_same():
  """Test with a list of two identical plugins with the same name. The first one should be erased."""
  from .plugins.mutation.snp_wrapper import SNP

  snp_params = {
          "model_id": "snp_test",
          "start_snps_frac": 0.1,
          "stop_snps_frac": 0.3,
          "phet": 0,
          "p": 0.01,
          "het_rng_seed": 3,
          "strand_rng_seed": 4,
          "poisson_rng_seed": 1,
          "base_sub_rng_seed": 2
  }
  snp_inputs = {}
  wrp_snp1 = SNP(snp_inputs, snp_params)
  wrp_snp1.execute()

  snp_params = {
          "model_id": "snp_test",
          "start_snps_frac": 0.7,
          "stop_snps_frac": 0.9,
          "phet": 1,
          "p": 0.01,
          "het_rng_seed": 3,
          "strand_rng_seed": 4,
          "poisson_rng_seed": 1,
          "base_sub_rng_seed": 2
  }
  snp_inputs = {}
  wrp_snp2 = SNP(snp_inputs, snp_params)
  wrp_snp2.execute()

  params = {'chromosome': '1'}
  inputs = {'ref': '/sbgenomics/test-data/porcine_circovirus.smalla',
            'plugins': [wrp_snp1.outputs.json_fragment, wrp_snp2.outputs.json_fragment]}
  wrp_m = Mutate(inputs, params)
  outputs = wrp_m.test()

  expected_line = '1	590	.	G	A	96	PASS	.	GT	0/1'
  with open(outputs.vcf[0], 'r') as f:
    assert f.readlines()[5].strip() == expected_line


def test_mutate_plugin_list_model_ids_different():
  """Test with a list of two identical plugins with the different model ids. Both the plugins
  should work as expected."""
  from .plugins.mutation.snp_wrapper import SNP
  from nose.tools import assert_equals

  snp_params = {
          "model_id": "snp_test",
          "start_snps_frac": 0.1,
          "stop_snps_frac": 0.3,
          "phet": 0,
          "p": 0.01,
          "het_rng_seed": 3,
          "strand_rng_seed": 4,
          "poisson_rng_seed": 1,
          "base_sub_rng_seed": 2
  }
  snp_inputs = {}
  wrp_snp1 = SNP(snp_inputs, snp_params)
  wrp_snp1.execute()

  snp_params = {
          "model_id": "snp_test2",
          "start_snps_frac": 0.7,
          "stop_snps_frac": 0.9,
          "phet": 1,
          "p": 0.01,
          "het_rng_seed": 3,
          "strand_rng_seed": 4,
          "poisson_rng_seed": 1,
          "base_sub_rng_seed": 2
  }
  snp_inputs = {}
  wrp_snp2 = SNP(snp_inputs, snp_params)
  wrp_snp2.execute()

  params = {'chromosome': '1'}
  inputs = {'ref': '/sbgenomics/test-data/porcine_circovirus.smalla',
            'plugins': [wrp_snp1.outputs.json_fragment, wrp_snp2.outputs.json_fragment]}
  wrp_m = Mutate(inputs, params)
  outputs = wrp_m.test()

  expected_lines = ['1	169	.	A	T	96	PASS	.	GT	1/1', '1	590	.	G	A	96	PASS	.	GT	0/1']
  with open(outputs.vcf[0], 'r') as f:
    lines = f.readlines()
    assert_equals(lines[5].strip(), expected_lines[0])
    assert_equals(lines[6].strip(), expected_lines[1])
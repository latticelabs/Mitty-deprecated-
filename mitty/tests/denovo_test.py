from nose.tools import assert_sequence_equal
from mitty.lib.variation import Variation, new_variation, HOMOZYGOUS, HET_10, HET_01

from mitty.lib.genome import FastaGenome
from mitty.lib.denovo import merge_variants_from_models, load_variant_model_list, main
from mitty.tests import example_fasta_genome


def add_variant_model_to_genome_test1():
  """Proposed data should collide and be ignored."""
  c1 = [v10, v11, v12] = [new_variation(12, 13, 'G', 'C', HET_01),
                          new_variation(13, 14, 'A', 'C', HET_10),
                          new_variation(19, 20, 'G', 'A', HET_10)]  # This should be placed with no problems
  g1 = {1: c1}
  vg = ({1: [new_variation(13, 14, 'A', 'T', HET_10)]} for _ in [1])  # This will collide with v11 and not be inserted
  g2 = merge_variants_from_models(g1, [vg])

  assert_sequence_equal(g1[1], g2[1])


def add_variant_model_to_genome_test2():
  """Variant generator adds new chromosomes."""
  g1 = {1: [new_variation(13, 14, 'A', 'C', HET_10)]}  # This should be placed with no problems
  vg = ({2: [new_variation(13, 14, 'A', 'T', HET_10)]} for _ in [1])
  g2 = merge_variants_from_models(g1, [vg])

  assert 1 in g2
  assert 2 in g2


def add_variant_model_to_genome_test3():
  """Variants from multiple models."""
  c1 = [v12, v13, v19] = [new_variation(12, 13, 'G', 'C', HOMOZYGOUS),
                          new_variation(13, 14, 'A', 'C', HOMOZYGOUS),
                          new_variation(19, 20, 'G', 'A', HOMOZYGOUS)]  # This should be placed with no problems
  g1 = {1: c1}
  v16 = new_variation(16, 17, 'T', 'C', HOMOZYGOUS)
  vg1 = ({1: [v12, v16]} for _ in [1])  # First data will collide
  v22 = new_variation(22, 23, 'T', 'C', HOMOZYGOUS)
  vg2 = ({1: [v22], 2: [v13, v19]} for _ in [1])  # This will all pass
  g2_correct = {
    1: [v12, v13, v16, v19, v22],
    2: [v13, v19]
  }
  g2 = merge_variants_from_models(g1, [vg1, vg2])

  assert_sequence_equal(g2[1], g2_correct[1])
  assert_sequence_equal(g2[2], g2_correct[2])


# This test uses the stock SNP which must exist for this test to pass
def load_variant_models_test():
  """Loading SNP data as a test."""
  param_json = [
        {
          "snp": {
             "phet": 0.5,
             "p": 0.01,
             "master_seed": 1
          }
        },
        {
          "snp": {
             "phet": 0.0,
             "p": 0.01,
             "master_seed": 2
          }
        }
      ]
  mdl = load_variant_model_list(param_json)
  assert hasattr(mdl[0]["model"], 'variant_generator')
  assert 'master_seed' in mdl[1]["params"]


# This test uses the stock SNP which must be functional for this test to pass as well as some files created during
# test setup
def main_test():
  """Integrating denovo and stock SNP"""
  import tempfile
  param_json = {
    "denovo_variant_models": [
      {
        "snp": {
           "chromosome": [1, 2],
           "phet": 0.5,
           "p": 0.0001
        }
      },
      {
        "snp": {
           "chromosome": [2],
           "phet": 0.0,
           "p": 0.01
        }
      }
    ]
  }
  _, vcf_file_fname = tempfile.mkstemp(suffix='.vcf.gz')
  ref = FastaGenome(seq_dir=example_fasta_genome)
  mdl = load_variant_model_list(param_json['denovo_variant_models'])
  g1 = main(ref, models=mdl, master_seed=1)
  assert type(g1) == dict
  assert type(g1[1][0]) == Variation
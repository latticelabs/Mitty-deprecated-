from mitty.denovo import *
from mitty.tests import *
from nose.tools import eq_, ok_


def add_variant_model_to_genome_test1():
  g1 = {1: [Variation(12, 13, 'G', 'C', HET_01),
            Variation(13, 14, 'A', 'C', HET_10),
            Variation(19, 20, 'G', 'A', HET_10)]}  # This should be placed with no problems

  def variant_generator():
    g2 = [{1: [Variation(13, 14, 'A', 'T', HET_10)]}]  # This will pass
    for g in g2:
      yield g

  correct_final_g = {1: [Variation(12, 13, 'G', 'C', HET_01),
                         Variation(13, 14, 'A', 'C', HET_10),
                         Variation(19, 20, 'G', 'A', HET_10)]}

  vg = variant_generator()
  g2 = merge_variants_with_genome(g1, vg)

  eq_(correct_final_g, g2)


def add_variant_model_to_genome_test1a():
  g1 = {1: [Variation(13, 14, 'A', 'C', HET_10)]}  # This should be placed with no problems

  def variant_generator():
    g2 = [{1: [Variation(13, 14, 'A', 'T', HET_10)]}]  # This will pass
    for g in g2:
      yield g

  correct_final_g = {1: [Variation(13, 14, 'A', 'C', HET_10)]}

  vg = variant_generator()
  g2 = merge_variants_with_genome(g1, vg)

  eq_(correct_final_g, g2)


def add_variant_model_to_genome_test2():
  g1 = {1: [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)],
        2: [Variation(7, 10, 'CAA', 'C', HET_01)]}  # This should be placed with no problems

  def variant_generator():
    g2 = [{1: [Variation(1, 2, 'C', 'CAA', HET_01)]},  # This will collide
          {2: [Variation(7, 8, 'G', 'T', HET_10)]},  # This will pass
          {2: [Variation(17, 18, 'G', 'T', HET_10)]}]  # This will pass
    for g in g2:
      yield g

  correct_final_g = {1: [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)],
                     2: [Variation(7, 10, 'CAA', 'C', HET_01),
                         Variation(7, 8, 'G', 'T', HET_10),
                         Variation(17, 18, 'G', 'T', HET_10)]}

  vg = variant_generator()
  g2 = merge_variants_with_genome(g1, vg)

  eq_(correct_final_g, g2)
  eq_(g1, {1: [Variation(1, 4, 'CAA', 'C', HOMOZYGOUS)], 2: [Variation(7, 10, 'CAA', 'C', HET_01)]})  # Must not change original


# This test uses the stock SNP which must exist for this test to pass
def load_variant_models_test():
  """Loading SNP variant as a test."""
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
  mdl = load_variant_models(param_json)
  assert hasattr(mdl[0]["model"], 'variant_generator')
  assert 'master_seed' in mdl[1]["params"]


# This test uses the stock SNP which must be functional for this test to pass as well as some files created during
# test setup
def main_test():
  """Integrating denovo and stock SNP"""
  import json, vcf, tempfile
  param_json = {
    "denovo_variant_models": [
      {
        "snp": {
           "chromosome": [1, 2],
           "phet": 0.5,
           "p": 0.0001,
           "master_seed": 1
        }
      },
      {
        "snp": {
           "chromosome": [2],
           "phet": 0.0,
           "p": 0.01,
           "master_seed": 2
        }
      }
    ]
  }
  _, vcf_file_fname = tempfile.mkstemp(suffix='.vcf.gz')
  ref = FastaGenome(seq_dir=example_fasta_genome)
  g1 = main(ref, vcf_file_name=vcf_file_fname, parameters=param_json, master_seed=1)
  assert os.path.exists(vcf_file_fname)

  vcf_rdr = vcf.Reader(filename=vcf_file_fname).fetch(chrom=1, start=0)
  v1 = vcf_rdr.next()
  eq_(v1.POS, g1[1][0].POS)